from flask import session, request, redirect, current_app, after_this_request, has_request_context, g
from server.auth.auth import AuthTypeClientBase, AuthTypeFactory
from server.common.errors import AuthenticationError, ConfigurationError
from urllib.parse import urlencode
from urllib.request import urlopen
import json

# It is not required to have authlib or jose.
# However, it is a configuration error to use this auth type if they are not installed.
missingimport = []
try:
    from authlib.integrations.flask_client import OAuth
except ModuleNotFoundError:
    missingimport.append("authlib")

try:
    from jose import jwt
except ModuleNotFoundError:
    missingimport.append("jose")


class AuthTypeOAuth(AuthTypeClientBase):
    """An authentication type for oauth2 logins."""

    CXG_ID_TOKEN = "id_token"

    def __init__(self, server_config):
        super().__init__()
        if missingimport:
            raise ConfigurationError(f"oauth requires these modules: {', '.join(missingimport)}")
        self.algorithms = ["RS256"]
        self.api_base_url = server_config.authentication__params_oauth__api_base_url
        self.client_id = server_config.authentication__params_oauth__client_id
        self.client_secret = server_config.authentication__params_oauth__client_secret
        self.callback_base_url = server_config.authentication__params_oauth__callback_base_url
        self.session_cookie = server_config.authentication__params_oauth__session_cookie
        self.cookie_params = server_config.authentication__params_oauth__cookie
        self._validate_cookie_params()

        # set the audience
        self.audience = self.client_id

        # load the jwks (JSON Web Key Set).
        # The JSON Web Key Set (JWKS) is a set of keys which contains the public keys used to verify
        # any JSON Web Token (JWT) issued by the authorization server and signed using the RS256
        try:
            jwksloc = f"{self.api_base_url}/.well-known/jwks.json"
            jwksurl = urlopen(jwksloc)
            self.jwks = json.loads(jwksurl.read())
        except Exception:
            raise ConfigurationError(f"error in oauth, api_url_base: {self.api_base_url}, cannot access {jwksloc}")

    def _validate_cookie_params(self):
        """check the cookie_params, and raise a ConfigurationError if there is something wrong"""
        if self.session_cookie:
            return

        if not isinstance(self.cookie_params, dict):
            raise ConfigurationError("either session_cookie or cookie must be set")
        valid_keys = {"key", "max_age", "expires", "path", "domain", "secure", "httponly", "samesite"}
        keys = set(self.cookie_params.keys())
        unknown = keys - valid_keys
        if unknown:
            raise ConfigurationError(f"unexpected key in cookie params: {', '.join(unknown)}")
        if "key" not in keys:
            raise ConfigurationError("must have a key (name) in the cookie params")

    def is_valid_authentication_type(self):
        return True

    def requires_client_login(self):
        return True

    def add_url_rules(self, app):
        app.add_url_rule("/login", "login", self.login, methods=["GET"])
        app.add_url_rule("/logout", "logout", self.logout, methods=["GET"])
        app.add_url_rule("/oauth2/callback", "callback", self.callback, methods=["GET"])

    def complete_setup(self, flask_app):
        self.oauth = OAuth(flask_app)
        if self.callback_base_url is None:
            # In this case, assume the server is running on the same host as the client,
            # and the oauth provider has been configured
            # with a callback that understands a localhost callback (e.g. A http://localhost:5005).
            server_config = flask_app.app_config.server_config
            self.callback_base_url = f"http://{server_config.app__host}:{server_config.app__port}"

        self.client = self.oauth.register(
            "oauth",
            client_id=self.client_id,
            client_secret=self.client_secret,
            api_base_url=self.api_base_url,
            access_token_url=f"{self.api_base_url}/oauth/token",
            authorize_url=f"{self.api_base_url}/authorize",
            client_kwargs={
                "scope" : "openid profile email",
            }
        )

    def is_user_authenticated(self):
        try:
            payload = self.get_jwt_payload()
            return payload is not None
        except AuthenticationError:
            return False

    def get_user_id(self):
        payload = self.get_jwt_payload()
        if payload and payload.get("sub"):
            return payload.get("sub")
        return None

    def get_user_name(self):
        payload = self.get_jwt_payload()
        if payload and payload.get("name"):
            return payload.get("name")
        return None

    def get_user_email(self):
        payload = self.get_jwt_payload()
        if payload and payload.get("email"):
            return payload.get("email")
        return None

    def update_response(self, response):
        response.cache_control.update(
            dict(public=True, max_age=0, no_store=True, no_cache=True, must_revalidate=True))

    def login(self):
        callbackurl = f'{self.callback_base_url}/oauth2/callback'
        return_path = request.args.get("dataset", "")
        return_to = f"{self.callback_base_url}/{return_path}"
        # save the return path in the session cookie, accessed in the callback function
        session["oauth_callback_redirect"] = return_to
        response = self.client.authorize_redirect(redirect_uri=callbackurl)
        self.update_response(response)
        return response

    def logout(self):
        if self.session_cookie:
            if self.CXG_ID_TOKEN in session:
                del session[self.CXG_ID_TOKEN]
        else:
            @after_this_request
            def remove_cookie(response):
                response.set_cookie(self.cookie_params["key"], "", expires=0)
                self.update_response(response)
                return response

        params = {'returnTo' : self.callback_base_url, 'client_id' : self.client_id}
        response = redirect(self.client.api_base_url + '/v2/logout?' + urlencode(params))
        self.update_response(response)
        return response

    def callback(self):
        token = self.client.authorize_access_token()
        id_token = token.get("id_token")
        oauth_callback_redirect = session.pop("oauth_callback_redirect", "/")
        resp = redirect(oauth_callback_redirect)

        if self.session_cookie:
            session[self.CXG_ID_TOKEN] = id_token
        else:
            args = self.cookie_params.copy()
            del args["key"]
            try:
                resp.set_cookie(
                    self.cookie_params["key"],
                    id_token,
                    **args)
                g.token = id_token
            except Exception as e:
                raise AuthenticationError(f"unable to set_cookie {self.cookie_params}") from e

        self.update_response(resp)
        return resp

    def get_login_url(self, data_adaptor):
        """Return the url for the login route"""
        if current_app.app_config.is_multi_dataset():
            return f"/login?dataset={data_adaptor.uri_path}"
        else:
            return "/login"

    def get_logout_url(self, data_adaptor):
        """Return the url for the logout route"""
        return "/logout"

    def get_token(self):
        """Function to return the token"""
        if "token" in g:
            return g.token
        if self.session_cookie:
            g.token = session.get(self.CXG_ID_TOKEN)
        else:
            g.token = request.cookies.get(self.cookie_params["key"])

        return g.token

    def get_jwt_payload(self):
        if not has_request_context():
            return None

        token = self.get_token()
        if token is None:
            return None

        unverified_header = jwt.get_unverified_header(token)
        rsa_key = {}
        for key in self.jwks['keys']:
            if key['kid'] == unverified_header['kid']:
                rsa_key = {
                    'kty': key['kty'],
                    'kid': key['kid'],
                    'use': key['use'],
                    'n': key['n'],
                    'e': key['e']
                }
        if rsa_key:
            try:
                payload = jwt.decode(
                    token,
                    rsa_key,
                    algorithms=self.algorithms,
                    audience=self.audience,
                    issuer=self.api_base_url + "/"
                )
                return payload

            except jwt.JWTError as e:
                raise AuthenticationError(f"invalid signature: {str(e)}")
            except jwt.ExpiredSignatureError as e:
                raise AuthenticationError(f"token expired: {str(e)}")
            except jwt.JWTClaimsError as e:
                raise AuthenticationError(f"invalid claims {str(e)}")

        raise AuthenticationError("Unable to find the appropriate key")


AuthTypeFactory.register("oauth", AuthTypeOAuth)
