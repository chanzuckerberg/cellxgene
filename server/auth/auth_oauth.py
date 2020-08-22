from flask import session, request, redirect, current_app, after_this_request, has_request_context, g
from server.auth.auth import AuthTypeClientBase, AuthTypeFactory
from server.common.errors import AuthenticationError, ConfigurationError
from urllib.parse import urlencode
import json
import requests
import base64

# It is not required to have authlib or jose.
# However, it is a configuration error to use this auth type if they are not installed.
missingimport = []
try:
    from authlib.integrations.flask_client import OAuth
except ModuleNotFoundError:
    missingimport.append("authlib")

try:
    from jose import jwt
    from jose.exceptions import ExpiredSignatureError, JWTError, JWTClaimsError
except ModuleNotFoundError:
    missingimport.append("jose")


class Tokens:
    """Simple class to represent the tokens that are saved/restored from the cookie"""

    def __init__(self, access_token, id_token, refresh_token, expires_at):
        self.access_token = access_token
        self.id_token = id_token
        self.refresh_token = refresh_token
        self.expires_at = expires_at
        if not (access_token and id_token and refresh_token and expires_at):
            raise KeyError(str(self.__dict__))


class AuthTypeOAuth(AuthTypeClientBase):
    """An authentication type for oauth2 logins."""

    CXG_TOKENS = "auth_tokens"

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
            jwksurl = requests.get(jwksloc)
            self.jwks = jwksurl.json()
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
            "auth0",
            client_id=self.client_id,
            client_secret=self.client_secret,
            api_base_url=self.api_base_url,
            refresh_token_url=f"{self.api_base_url}/oauth/token",
            access_token_url=f"{self.api_base_url}/oauth/token",
            authorize_url=f"{self.api_base_url}/authorize",
            client_kwargs={"scope": "openid profile email offline_access"},
        )

    def is_user_authenticated(self):
        payload = self.get_userinfo()
        return payload is not None

    def get_user_id(self):
        payload = self.get_userinfo()
        if payload and payload.get("sub"):
            return payload.get("sub")
        return None

    def get_user_name(self):
        payload = self.get_userinfo()
        if payload and payload.get("name"):
            return payload.get("name")
        return None

    def get_user_email(self):
        payload = self.get_userinfo()
        if payload and payload.get("email"):
            return payload.get("email")
        return None

    def update_response(self, response):
        response.cache_control.update(dict(public=True, max_age=0, no_store=True, no_cache=True, must_revalidate=True))

    def login(self):
        callbackurl = f"{self.callback_base_url}/oauth2/callback"
        return_path = request.args.get("dataset", "")
        return_to = f"{self.callback_base_url}/{return_path}"
        # save the return path in the session cookie, accessed in the callback function
        session["oauth_callback_redirect"] = return_to
        response = self.client.authorize_redirect(redirect_uri=callbackurl)
        self.update_response(response)
        return response

    def logout(self):
        self.remove_tokens()
        params = {"returnTo": self.callback_base_url, "client_id": self.client_id}
        response = redirect(self.client.api_base_url + "/v2/logout?" + urlencode(params))
        self.update_response(response)
        return response

    def callback(self):
        data = self.client.authorize_access_token()
        tokens = Tokens(
            access_token=data.get("access_token"),
            id_token=data.get("id_token"),
            refresh_token=data.get("refresh_token"),
            expires_at=data.get("expires_at"),
        )
        self.save_tokens(tokens)
        oauth_callback_redirect = session.pop("oauth_callback_redirect", "/")
        response = redirect(oauth_callback_redirect)
        self.update_response(response)
        return response

    def get_tokens(self):
        """Extract the tokens from the cookie, and store them in the flask global context"""
        if "tokens" in g:
            return g.tokens

        try:
            if self.session_cookie:
                tokensdict = session.get(self.CXG_TOKENS)
                if tokensdict:
                    g.tokens = Tokens(**tokensdict)
                else:
                    return None
            else:
                value = request.cookies.get(self.cookie_params["key"])
                value = base64.b64decode(value)
                try:
                    tokensdict = json.loads(value)
                    g.tokens = Tokens(**tokensdict)
                except (TypeError, KeyError, json.decoder.JSONDecodeError):
                    g.pop("tokens", None)
                    return None

        except (TypeError, KeyError):
            g.pop("tokens", None)
            return None

        return g.tokens

    def save_tokens(self, tokens):
        g.tokens = tokens
        if self.session_cookie:
            session[self.CXG_TOKENS] = tokens.__dict__
        else:

            @after_this_request
            def set_cookie(response):
                args = self.cookie_params.copy()
                value = base64.b64encode(json.dumps(tokens.__dict__).encode("utf-8"))
                del args["key"]
                try:
                    response.set_cookie(self.cookie_params["key"], value, **args)
                except Exception as e:
                    raise AuthenticationError(f"unable to set_cookie {self.cookie_params}") from e
                return response

    def remove_tokens(self):
        g.pop("tokens", None)
        if self.session_cookie:
            if self.CXG_TOKENS in session:
                del session[self.CXG_TOKENS]
        else:

            @after_this_request
            def remove_cookie(response):
                response.set_cookie(self.cookie_params["key"], "", expires=0)
                self.update_response(response)
                return response

    def get_login_url(self, data_adaptor):
        """Return the url for the login route"""
        if current_app.app_config.is_multi_dataset():
            return f"/login?dataset={data_adaptor.uri_path}/"
        else:
            return "/login"

    def get_logout_url(self, data_adaptor):
        """Return the url for the logout route"""
        return "/logout"

    def check_jwt_payload(self, id_token):
        try:
            unverified_header = jwt.get_unverified_header(id_token)
        except JWTError:
            return None

        rsa_key = {}
        for key in self.jwks["keys"]:
            if key["kid"] == unverified_header["kid"]:
                rsa_key = {
                    "kty": key["kty"],
                    "kid": key["kid"],
                    "use": key["use"],
                    "n": key.get("n"),
                    "e": key.get("e"),
                }
        if rsa_key:
            options = {}
            if not rsa_key["n"] or not rsa_key["e"]:
                # this is a mock auth server, do not validate
                options = {"verify_signature": False, "verify_iss": False}
            try:
                payload = jwt.decode(
                    id_token,
                    rsa_key,
                    algorithms=self.algorithms,
                    audience=self.audience,
                    issuer=self.api_base_url + "/",
                    options=options,
                )
                return payload

            except ExpiredSignatureError:
                # This exception is handled in get_userinfo
                raise
            except JWTClaimsError as e:
                raise AuthenticationError(f"invalid claims {str(e)}") from e
            except JWTError as e:
                raise AuthenticationError(f"invalid signature: {str(e)}") from e

        raise AuthenticationError("Unable to find the appropriate key")

    def get_userinfo(self):
        if not has_request_context():
            return None

        # check if the userinfo has been retrieved already in this request
        if "userinfo" in g:
            return g.get("userinfo")

        # if there is no id_token, return None (user is not authenticated)
        tokens = self.get_tokens()
        if tokens is None or tokens.id_token is None:
            return None

        try:
            # check the jwt payload.  This raises an AuthenticationError if the token is not valid.
            # It the token has expired, we attempt to refresh the token
            g.userinfo = self.check_jwt_payload(tokens.id_token)
            return g.userinfo

        except ExpiredSignatureError:
            tokens = self.refresh_expired_token(tokens.refresh_token)
            if tokens is None or tokens.id_token is None:
                return None
            else:
                try:
                    g.userinfo = self.check_jwt_payload(tokens.id_token)
                    return g.userinfo
                except JWTError as e:
                    raise AuthenticationError(f"error during token refresh: {str(e)}") from e

        except AuthenticationError:
            self.remove_tokens()
            raise

    def refresh_expired_token(self, refresh_token):
        params = {
            "grant_type": "refresh_token",
            "client_id": self.client_id,
            "refresh_token": refresh_token,
            "client_secret": self.client_secret,
        }
        headers = {"content-type": "application/x-www-form-urlencoded"}
        request = requests.post(f"{self.api_base_url}/oauth/token", urlencode(params), headers=headers)
        if request.status_code != 200:
            # unable to refresh the token, log the user out
            self.remove_tokens()
            return None
        data = request.json()
        tokens = Tokens(
            access_token=data.get("access_token"),
            id_token=data.get("id_token"),
            refresh_token=data.get("refresh_token", refresh_token),
            expires_at=data.get("expires_at"),
        )
        self.save_tokens(tokens)
        return tokens


AuthTypeFactory.register("oauth", AuthTypeOAuth)
