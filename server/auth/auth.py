from abc import ABC, abstractmethod


class AuthTypeBase(ABC):
    """Base type for all authentication types."""

    def __init__(self):
        super().__init__()

    @abstractmethod
    def set_params(self, params):
        """Set the parameters from app config.  raise ConfigurationError if any params are invalid"""
        pass

    @abstractmethod
    def is_valid(self):
        """Return True if the auth type can return user info (AuthTypeNone is the only one that cannot)"""
        pass

    def requires_client_login(self):
        """Return True if the user needs to login from the client (e.g. Login button is shown)"""
        return False

    @abstractmethod
    def is_authenticated(self):
        """Return True if the user is authenticated"""
        pass

    @abstractmethod
    def get_userid(self):
        """Return the id for this user (string)"""
        pass

    @abstractmethod
    def get_username(self):
        """Return the name of the user (string)"""
        pass


class AuthTypeClientBase(AuthTypeBase):
    """Base type for all authentication types that require the client to login"""

    def __init__(self):
        super().__init__()

    def requires_client_login(self):
        return True

    @abstractmethod
    def add_url_rules(self, selfapp):
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
        assert(issubclass(auth_type, AuthTypeBase))
        AuthTypeFactory.auth_types[name] = auth_type

    @staticmethod
    def create(name):
        auth_type = AuthTypeFactory.auth_types.get(name)
        if auth_type is None:
            return None
        return auth_type()
