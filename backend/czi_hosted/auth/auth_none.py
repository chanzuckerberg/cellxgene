from backend.czi_hosted.auth.auth import AuthTypeBase, AuthTypeFactory


class AuthTypeNone(AuthTypeBase):
    def __init__(self, app_config):
        super().__init__()

    def is_valid_authentication_type(self):
        return False

    def complete_setup(self, app):
        pass

    def is_user_authenticated(self):
        return True

    def get_user_id(self):
        return None

    def get_user_name(self):
        return None

    def get_user_email(self):
        return None


AuthTypeFactory.register(None, AuthTypeNone)
