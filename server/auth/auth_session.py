from server.auth.auth import AuthTypeBase, AuthTypeFactory
from flask import session
from uuid import uuid4


class AuthTypeSession(AuthTypeBase):
    """Session based authentication.  The user is always logged.  The user id is a random number
    associated with the session.  This is a good choice for desktop servers."""

    # key in the session token for userid
    CXGUID = "cxguid"

    def __init__(self, app_config):
        super().__init__()

    def is_valid(self):
        return True

    def complete_setup(self, app):
        pass

    def is_authenticated(self):
        # always authenticated
        return True

    def get_userid(self):
        if self.CXGUID not in session:
            session[self.CXGUID] = uuid4().hex
            session.permanent = True
        return session[self.CXGUID]

    def get_username(self):
        return "anonymous"


AuthTypeFactory.register("session", AuthTypeSession)
