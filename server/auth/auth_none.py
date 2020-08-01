from server.auth.auth import AuthTypeBase, AuthTypeFactory
from server.common.errors import ConfigurationError


class AuthTypeNone(AuthTypeBase):
    def __init__(self):
        super().__init__()

    def is_valid(self):
        return False

    def set_params(self, params):
        if params:
            raise ConfigurationError("not expecting authentication parameters")

    def is_authenticated(self):
        return True

    def get_userid(self):
        return None

    def get_username(self):
        return None


AuthTypeFactory.register(None, AuthTypeNone)
