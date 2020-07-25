from server.auth.auth import AuthTypeBase, AuthTypeFactory


class AuthTypeNone(AuthTypeBase):

    def __init__(self, app_config):
        super().__init__()

    def is_valid(self):
        return False

    def complete_setup(self, app):
        pass

    def is_authenticated(self):
        return True

    def get_userid(self):
        return None

    def get_username(self):
        return None


AuthTypeFactory.register(None, AuthTypeNone)
