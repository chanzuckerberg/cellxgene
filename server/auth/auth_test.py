from server.auth.auth import AuthTypeClientBase, AuthTypeFactory
from flask import session, request, redirect, current_app


class AuthTypeTest(AuthTypeClientBase):
    """An authentication type for testing client based logins.  When the login route is accessed
    the user is automatically logged in with a default or configured username"""

    # key in session token with userid and username
    CXGUID = "cxguid_test"
    CXGUNAME = "cxguname_test"
    CXGUEMAIL = "cxguemail_test"
    CXGUPICTURE = "cxgupicture_test"

    def __init__(self, app_config):
        super().__init__()
        self.user_name = "test_account"
        self.user_id = "id0001"
        self.user_email = "test_account@test.com"
        self.user_picture = None

    def is_valid_authentication_type(self):
        return True

    def requires_client_login(self):
        return True

    def add_url_rules(self, app):
        app.add_url_rule("/login", "login", self.login, methods=["GET"])
        app.add_url_rule("/logout", "logout", self.logout, methods=["GET"])

    def complete_setup(self, app):
        pass

    def is_user_authenticated(self):
        return self.CXGUID in session

    def get_user_id(self):
        return session.get(self.CXGUID)

    def get_user_name(self):
        return session.get(self.CXGUNAME)

    def get_user_email(self):
        return session.get(self.CXGUEMAIL)

    def get_user_picture(self):
        return session.get(self.CXGUPICTURE)

    def login(self):
        args = request.args
        return_to = args.get("dataset", "/")
        session[self.CXGUID] = args.get("userid", self.user_id)
        session[self.CXGUNAME] = args.get("username", self.user_name)
        session[self.CXGUEMAIL] = args.get("email", self.user_email)
        session[self.CXGUPICTURE] = args.get("picture", self.user_picture)
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
