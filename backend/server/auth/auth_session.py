from flask import session
from uuid import uuid4


class Session:
    """ Session based tracking """

    # key in the session token for userid
    CXGUID = "cxguid"

    def get_user_id(self):
        if self.CXGUID not in session:
            session[self.CXGUID] = uuid4().hex
            session.permanent = True
        return session[self.CXGUID]
