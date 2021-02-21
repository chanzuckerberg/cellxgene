from abc import ABC, abstractmethod


class AuthTypeBase(ABC):
    """Base type for all authentication types."""

    def __init__(self):
        super().__init__()

    @abstractmethod
    def is_valid_authentication_type(self):
        """Return True if the auth type is valid, e.g. it can return userinfo and username.
          (AuthTypeNone is the only one type that returns False)"""
        pass

    def requires_client_login(self):
        """Return True if the user needs to login from the client (e.g. Login button is shown)"""
        return False

    @abstractmethod
    def complete_setup(self, app):
        """complete any setup that may be needed by this auth type.  The Flask app is passed in.
        This is the last auth function called before the server starts to run."""
        pass

    @abstractmethod
    def is_user_authenticated(self):
        """Return True if the user is authenticated"""
        pass

    @abstractmethod
    def get_user_id(self):
        """Return the id for this user (string)"""
        pass

    @abstractmethod
    def get_user_name(self):
        """Return the name of the user (string)"""
        pass

    @abstractmethod
    def get_user_email(self):
        """Return the name of the user (string)"""
        pass

    def get_user_picture(self):
        """Return the location to the user's picture"""
        return None


class AuthTypeClientBase(AuthTypeBase):
    """Base type for all authentication types that require the client to login"""

    def __init__(self):
        super().__init__()

    def requires_client_login(self):
        return True

    @abstractmethod
    def add_url_rules(self, app):
        """Add url rules to the app (like /login, /logout, etc)"""
        pass

    @abstractmethod
    def get_login_url(self, data_adaptor):
        """Return the url for the login route"""
        pass

    @abstractmethod
    def get_logout_url(self, data_adaptor):
        """Return the url for the logout route"""
        pass


class AuthTypeFactory:
    """Factory class to create an authentication type"""

    auth_types = {}

    @staticmethod
    def register(name, auth_type):
        assert issubclass(auth_type, AuthTypeBase)
        AuthTypeFactory.auth_types[name] = auth_type

    @staticmethod
    def create(name, app_config):
        auth_type = AuthTypeFactory.auth_types.get(name)
        if auth_type is None:
            return None
        return auth_type(app_config)
