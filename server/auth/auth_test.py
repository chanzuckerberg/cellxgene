from server.auth.auth import AuthTypeClientBase, AuthTypeFactory
from flask import session, request, redirect, current_app


class AuthTypeTest(AuthTypeClientBase):
    """An authentication type for testing client based logins.  When the login route is accessed
    the user is automatically logged in with a default or configured username"""

    # key in session token with userid and username
    CXGUID = "cxguid_test"
    CXGUNAME = "cxguname_test"

    def __init__(self):
        super().__init__()
        self.username = "test_account"
        self.userid = "id0001"

    def is_valid(self):
        return True

    def requires_client_login(self):
        return True

    def add_url_rules(self, app):
        app.add_url_rule("/login", "login", self.login, methods=["GET"])
        app.add_url_rule("/logout", "logout", self.logout, methods=["GET"])

    def set_params(self, params):
        if params:
            self.username = params.get("username", self.username)
            self.userid = params.get("userid", self.userid)

    def is_authenticated(self):
        return self.CXGUID in session

    def get_userid(self):
        return session.get(self.CXGUID)

    def get_username(self):
        return session.get(self.CXGUNAME)

    def login(self):
        args = request.args
        return_to = args.get("dataset", "/")
        session[self.CXGUID] = args.get("userid", self.userid)
        session[self.CXGUNAME] = args.get("username", self.username)
        return redirect(return_to)

    def logout(self):
        session.clear()
        return_to = request.args.get("dataset", "/")
        return redirect(return_to)

    def get_login_url(self, data_adaptor):
        """Return the url for the login route"""
        if current_app.app_config.is_multi_dataset():
            return f"/login?dataset={data_adaptor.uri_path}"
        else:
            return "/login"

    def get_logout_url(self, data_adaptor):
        """Return the url for the logout route"""
        if current_app.app_config.is_multi_dataset():
            return f"/logout?dataset={data_adaptor.uri_path}"
        else:
            return "/logout"


AuthTypeFactory.register("test", AuthTypeTest)
