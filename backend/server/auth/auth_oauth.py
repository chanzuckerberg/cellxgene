from flask import session, request, redirect, current_app, after_this_request, has_request_context, g
from backend.server.auth.auth import AuthTypeClientBase, AuthTypeFactory
from backend.server.common.errors import AuthenticationError, ConfigurationError
from urllib.parse import urlencode, urlparse
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

    def __init__(self, access_token, id_token, refresh_token, expires_at, **kwargs):
        self.access_token = access_token
        self.id_token = id_token
        self.refresh_token = refresh_token
        self.expires_at = expires_at

        # expires_at may be None after a token refresh, and so it is not checked here
        if not (access_token and id_token and refresh_token):
            raise KeyError(str(self.__dict__))


class AuthTypeOAuth(AuthTypeClientBase):
    """An authentication type for oauth2 logins."""

    CXG_TOKENS = "auth_tokens"

    def __init__(self, server_config):
        super().__init__()
        if missingimport:
            raise ConfigurationError(f"oauth requires these modules: {', '.join(missingimport)}")
        self.algorithms = ["RS256"]
        self.oauth_api_base_url = server_config.authentication__params_oauth__oauth_api_base_url
        self.client_id = server_config.authentication__params_oauth__client_id
        self.client_secret = server_config.authentication__params_oauth__client_secret
        self.session_cookie = server_config.authentication__params_oauth__session_cookie
        self.cookie_params = server_config.authentication__params_oauth__cookie
        self.jwt_decode_options = server_config.authentication__params_oauth__jwt_decode_options

        self._validate_cookie_params()
        self._validate_jwt_decode_options()

        self.api_base_url = server_config.get_api_base_url()
        self.web_base_url = server_config.get_web_base_url()
        if self.api_base_url is None:
            raise ConfigurationError("oauth requires the app__api_base_url to be set")

        # set the audience
        self.audience = self.client_id

        # load the jwks (JSON Web Key Set).
        # The JSON Web Key Set (JWKS) is a set of keys which contains the public keys used to verify
        # any JSON Web Token (JWT) issued by the authorization server and signed using the RS256
        try:
            jwksloc = f"{self.oauth_api_base_url}/.well-known/jwks.json"
            jwksurl = requests.get(jwksloc)
            self.jwks = jwksurl.json()
        except Exception:
            raise ConfigurationError(
                f"error in oauth, api_url_base: {self.oauth_api_base_url}, cannot access {jwksloc}"
            )

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

    def _validate_jwt_decode_options(self):
        """check the jwt_decode_options, and raise a ConfigurationError if there is something wrong"""
        if self.jwt_decode_options is None:
            self.jwt_decode_options = {}
            return

        valid_keys = {
            "verify_signature",
            "verify_aud",
            "verify_iat",
            "verify_exp",
            "verify_nbf",
            "verify_iss",
            "verify_sub",
            "verify_jti",
            "verify_at_hash",
            "leeway",
        }
        keys = set(self.jwt_decode_options.keys())
        unknown = keys - valid_keys
        if unknown:
            raise ConfigurationError(f"unexpected key in jwt_decode_options: {', '.join(unknown)}")

    def is_valid_authentication_type(self):
        return True

    def requires_client_login(self):
        return True

    def add_url_rules(self, app):
        parse = urlparse(self.api_base_url)
        app.add_url_rule(f"{parse.path}/login", "login", self.login, methods=["GET"])
        app.add_url_rule(f"{parse.path}/logout", "logout", self.logout, methods=["GET"])
        app.add_url_rule(f"{parse.path}/logout_redirect", "logout_redirect", self.logout_redirect, methods=["GET"])
        app.add_url_rule(f"{parse.path}/oauth2/callback", "callback", self.callback, methods=["GET"])

    def complete_setup(self, flask_app):
        self.oauth = OAuth(flask_app)

        self.client = self.oauth.register(
            "auth0",
            client_id=self.client_id,
            client_secret=self.client_secret,
            api_base_url=self.oauth_api_base_url,
            refresh_token_url=f"{self.oauth_api_base_url}/oauth/token",
            access_token_url=f"{self.oauth_api_base_url}/oauth/token",
            authorize_url=f"{self.oauth_api_base_url}/authorize",
            client_kwargs={"scope": "openid profile email offline_access"},
        )

    def is_user_authenticated(self):
        payload = self.get_userinfo()
        return payload is not None

    def get_user_id(self):
        payload = self.get_userinfo()
        return payload.get("sub") if payload else None

    def get_user_name(self):
        payload = self.get_userinfo()
        return payload.get("name") if payload else None

    def get_user_email(self):
        payload = self.get_userinfo()
        return payload.get("email") if payload else None

    def get_user_picture(self):
        payload = self.get_userinfo()
        return payload.get("picture") if payload else None

    def update_response(self, response):
        response.cache_control.update(dict(public=True, max_age=0, no_store=True, no_cache=True, must_revalidate=True))

    def login(self):
        callbackurl = f"{self.api_base_url}/oauth2/callback"
        return_path = request.args.get("dataset", "")
        return_to = f"{self.web_base_url}/{return_path}"
        # save the return path in the session cookie, accessed in the callback function
        session["oauth_callback_redirect"] = return_to
        response = self.client.authorize_redirect(redirect_uri=callbackurl)
        self.update_response(response)
        return response

    def logout(self):
        """
        We would like for the user to remain on the same dataset after logout.  oauth requires that
        the redirect `returnTo` path be whitelisted by the oauth server, therefore a level of
        indirection is used.  We first redirect to a single path "logout_redirect", and logout_redirect
        will redirect the user's browser back to the current page.
        """
        self.remove_tokens()
        redirect_path = request.args.get("dataset", "")
        redirect_to = f"{self.web_base_url}/{redirect_path}"
        session["oauth_logout_redirect"] = redirect_to

        return_to = f"{self.api_base_url}/logout_redirect"
        params = {"returnTo": return_to, "client_id": self.client_id}
        response = redirect(self.client.api_base_url + "/v2/logout?" + urlencode(params))
        self.update_response(response)
        return response

    def logout_redirect(self):
        oauth_logout_redirect = session.pop("oauth_logout_redirect", "/")
        response = redirect(oauth_logout_redirect)
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
                value = session.get(self.CXG_TOKENS)
                if value:
                    g.tokens = Tokens(**value)
                else:
                    return None
            else:
                value = request.cookies.get(self.cookie_params["key"])
                if value is None:
                    return None
                value = base64.b64decode(value)
                value = json.loads(value)
                g.tokens = Tokens(**value)

        except Exception:
            # there are many types of exceptions that can be raise in the above section.
            # It is impractical to list all the exceptions here, since that would be brittle.
            # If an exception occurs, then return None, meaning that no token could be retrieved.
            current_app.logger.warning(f"auth cookie is in the wrong format: {str(value)}")
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
        if data_adaptor and current_app.app_config.is_multi_dataset():
            return f"{self.api_base_url}/login?dataset={data_adaptor.uri_path}/"
        else:
            return f"{self.api_base_url}/login"

    def get_logout_url(self, data_adaptor):
        """Return the url for the logout route"""
        if data_adaptor and current_app.app_config.is_multi_dataset():
            return f"{self.api_base_url}/logout?dataset={data_adaptor.uri_path}/"
        else:
            return f"{self.api_base_url}/logout"

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
            try:
                payload = jwt.decode(
                    id_token,
                    rsa_key,
                    algorithms=self.algorithms,
                    audience=self.audience,
                    issuer=self.oauth_api_base_url + "/",
                    options=self.jwt_decode_options,
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
        request = requests.post(f"{self.oauth_api_base_url}/oauth/token", urlencode(params), headers=headers)
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
