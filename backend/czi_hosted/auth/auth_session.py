from flask import session
from uuid import uuid4

from backend.czi_hosted.auth.auth import AuthTypeBase, AuthTypeFactory


class AuthTypeSession(AuthTypeBase):
    """Session based authentication.  The user is always logged.  The user id is a random number
    associated with the session.  This is a good choice for desktop servers."""

    # key in the session token for userid
    CXGUID = "cxguid"

    def __init__(self, app_config):
        super().__init__()

    def is_valid_authentication_type(self):
        return True

    def complete_setup(self, app):
        pass

    def is_user_authenticated(self):
        # always authenticated
        return True

    def get_user_id(self):
        if self.CXGUID not in session:
            session[self.CXGUID] = uuid4().hex
            session.permanent = True
        return session[self.CXGUID]

    def get_user_name(self):
        return "anonymous"

    def get_user_email(self):
        return None


AuthTypeFactory.register("session", AuthTypeSession)
