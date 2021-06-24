import json
import unittest

from backend.czi_hosted.common.config.app_config import AppConfig
from backend.test import FIXTURES_ROOT
from backend.test.test_czi_hosted.unit import BaseTest


class AuthTest(BaseTest):
    def setUp(self):
        self.dataset_dataroot = FIXTURES_ROOT

    def test_auth_none(self):
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(authentication__type=None, multi_dataset__dataroot=self.dataset_dataroot)
        app_config.update_default_dataset_config(user_annotations__enable=False)

        app_config.complete_config()
        server= self.create_app(app_config)
        server.testing = True
        session = server.test_client()
        config = json.loads(session.get(f"{self.TEST_URL_BASE}config").data)
        userinfo = json.loads(session.get(f"{self.TEST_URL_BASE}userinfo").data)
        self.assertNotIn("authentication", config["config"])
        self.assertIsNone(userinfo)

    def test_auth_session(self):
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(authentication__type="session", multi_dataset__dataroot=self.dataset_dataroot)
        app_config.update_default_dataset_config(user_annotations__enable=True)
        app_config.complete_config()

        server = self.create_app(app_config)
        server.auth.is_user_authenticated = lambda: True
        server.testing = True
        session = server.test_client()
        config = json.loads(session.get(f"{self.TEST_URL_BASE}config").data)
        userinfo = json.loads(session.get(f"{self.TEST_URL_BASE}userinfo").data)

        self.assertFalse(config["config"]["authentication"]["requires_client_login"])
        self.assertTrue(userinfo["userinfo"]["is_authenticated"])
        self.assertEqual(userinfo["userinfo"]["username"], "anonymous")

    def test_auth_test(self):
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(authentication__type="test")
        app_config.update_server_config(authentication__insecure_test_environment=True)
        app_config.update_server_config(
            multi_dataset__dataroot=dict(
                a1=dict(dataroot=self.dataset_dataroot, base_url="auth"),
                a2=dict(dataroot=self.dataset_dataroot, base_url="no-auth"),
            )
        )

        # specialize the configs
        app_config.add_dataroot_config("a1", app__authentication_enable=True, user_annotations__enable=True)
        app_config.add_dataroot_config("a2", app__authentication_enable=False, user_annotations__enable=False)

        app_config.complete_config()

        server=self.create_app(app_config)
        server.testing = True
        session = server.test_client()

        # auth datasets
        config = json.loads(session.get("/auth/pbmc3k.cxg/api/v0.2/config").data)
        userinfo = json.loads(session.get(f"/auth/pbmc3k.cxg/api/v0.2/userinfo").data)

        self.assertFalse(userinfo["userinfo"]["is_authenticated"])
        self.assertIsNone(userinfo["userinfo"]["username"])
        self.assertTrue(config["config"]["authentication"]["requires_client_login"])
        self.assertTrue(config["config"]["parameters"]["annotations"])

        login_uri = config["config"]["authentication"]["login"]
        logout_uri = config["config"]["authentication"]["logout"]

        self.assertEqual(login_uri, "/login?dataset=auth/pbmc3k.cxg")
        self.assertEqual(logout_uri, "/logout?dataset=auth/pbmc3k.cxg")

        response = session.get(login_uri)
        # check that the login redirect worked

        self.assertEqual(response.status_code, 302)
        self.assertEqual(response.headers['Location'], 'http://localhost/auth/pbmc3k.cxg')
        config = json.loads(session.get("/auth/pbmc3k.cxg/api/v0.2/config").data)
        userinfo = json.loads(session.get("/auth/pbmc3k.cxg/api/v0.2/userinfo").data)

        self.assertTrue(userinfo["userinfo"]["is_authenticated"])
        self.assertEqual(userinfo["userinfo"]["username"], "test_account")
        self.assertEqual(userinfo["userinfo"]["picture"], None)
        self.assertTrue(config["config"]["parameters"]["annotations"])

        response = session.get(logout_uri)
        # check that the logout redirect worked

        self.assertEqual(response.status_code, 302)
        config = json.loads(session.get("/auth/pbmc3k.cxg/api/v0.2/config").data)
        userinfo = json.loads(session.get("/auth/pbmc3k.cxg/api/v0.2/userinfo").data)
        self.assertFalse(userinfo["userinfo"]["is_authenticated"])
        self.assertIsNone(userinfo["userinfo"]["username"])
        self.assertTrue(config["config"]["parameters"]["annotations"])

        # no-auth datasets
        config = json.loads(session.get("/no-auth/pbmc3k.cxg/api/v0.2/config").data)
        userinfo = json.loads(session.get("/no-auth/pbmc3k.cxg/api/v0.2/userinfo").data)
        self.assertIsNone(userinfo)
        self.assertFalse(config["config"]["parameters"]["annotations"])

        # login with a picture
        session.get(f"{login_uri}&picture=myimage.png")
        userinfo = json.loads(session.get("/auth/pbmc3k.cxg/api/v0.2/userinfo").data)
        self.assertTrue(userinfo["userinfo"]["is_authenticated"])
        self.assertEqual(userinfo["userinfo"]["picture"], "myimage.png")

    def test_auth_test_single(self):
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(
            authentication__type="test", single_dataset__datapath=f"{self.dataset_dataroot}/pbmc3k.cxg"
        )
        app_config.update_server_config(authentication__insecure_test_environment=True)

        app_config.complete_config()

        server = self.create_app(app_config)
        server.testing = True
        session = server.test_client()

        config = json.loads(session.get("/api/v0.2/config").data)
        userinfo = json.loads(session.get("/api/v0.2/userinfo").data)
        self.assertFalse(userinfo["userinfo"]["is_authenticated"])
        self.assertIsNone(userinfo["userinfo"]["username"])
        self.assertTrue(config["config"]["authentication"]["requires_client_login"])
        self.assertTrue(config["config"]["parameters"]["annotations"])

        login_uri = config["config"]["authentication"]["login"]
        logout_uri = config["config"]["authentication"]["logout"]

        self.assertEqual(login_uri, "/login")
        self.assertEqual(logout_uri, "/logout")


        # check that the login redirect worked
        with server.test_client() as session:
            response = session.get(login_uri)
            self.assertEqual(response.status_code, 302)
            self.assertEqual(response.headers['Location'], "http://localhost/")

        config = json.loads(session.get("api/v0.2/config").data)
        userinfo = json.loads(session.get("/api/v0.2/userinfo").data)
        self.assertTrue(userinfo["userinfo"]["is_authenticated"])
        self.assertEqual(userinfo["userinfo"]["username"], "test_account")
        self.assertTrue(config["config"]["parameters"]["annotations"])

        response = session.get(logout_uri)
        # check that the logout redirect worked

        self.assertEqual(response.status_code, 302)
        self.assertEqual(response.headers['Location'], "http://localhost/")
        config = json.loads(session.get("/api/v0.2/config").data)

        userinfo = json.loads(session.get("/api/v0.2/userinfo").data)
        self.assertFalse(userinfo["userinfo"]["is_authenticated"])
        self.assertIsNone(userinfo["userinfo"]["username"])
        self.assertTrue(config["config"]["parameters"]["annotations"])
