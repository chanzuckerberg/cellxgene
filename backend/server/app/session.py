from uuid import uuid4

from flask.sessions import SessionMixin

CXGUID = "cxguid"


def get_user_id(session: SessionMixin) -> str:
    """ Gets a session-persistent user id. Creates one in the Flask session if non-extant """
    if CXGUID not in session:
        session[CXGUID] = uuid4().hex
        session.permanent = True
    return session[CXGUID]
