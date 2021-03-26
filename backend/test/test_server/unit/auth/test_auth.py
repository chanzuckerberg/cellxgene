import unittest

import requests

from backend.server.common.config.app_config import AppConfig
from backend.test.test_server.unit import test_server
from backend.test import H5AD_FIXTURE


class AuthTest(unittest.TestCase):
    def setUp(self):
        self.dataset_datapath = H5AD_FIXTURE

    def test_auth_none(self):
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(authentication__type=None, single_dataset__datapath=self.dataset_datapath)
        app_config.update_dataset_config(user_annotations__enable=False, user_annotations__gene_sets__readonly=True)

        app_config.complete_config()

        with test_server(app_config=app_config) as server:
            session = requests.Session()
            config = session.get(f"{server}/api/v0.2/config").json()
            userinfo = session.get(f"{server}/api/v0.2/userinfo").json()
            self.assertNotIn("authentication", config["config"])
            self.assertIsNone(userinfo)

    def test_auth_session(self):
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(authentication__type="session", single_dataset__datapath=self.dataset_datapath)
        app_config.update_dataset_config(user_annotations__enable=True)
        app_config.complete_config()
        with test_server(app_config=app_config) as server:
            session = requests.Session()
            config = session.get(f"{server}/api/v0.2/config").json()
            userinfo = session.get(f"{server}/api/v0.2/userinfo").json()

            self.assertFalse(config["config"]["authentication"]["requires_client_login"])
            self.assertTrue(userinfo["userinfo"]["is_authenticated"])
            self.assertEqual(userinfo["userinfo"]["username"], "anonymous")

    def test_auth_test_single(self):
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(
            authentication__type="test",
            single_dataset__datapath=self.dataset_datapath,
            authentication__insecure_test_environment=True,
        )

        app_config.complete_config()

        with test_server(app_config=app_config) as server:
            session = requests.Session()
            config = session.get(f"{server}/api/v0.2/config").json()
            userinfo = session.get(f"{server}/api/v0.2/userinfo").json()
            self.assertFalse(userinfo["userinfo"]["is_authenticated"])
            self.assertIsNone(userinfo["userinfo"]["username"])
            self.assertTrue(config["config"]["authentication"]["requires_client_login"])
            self.assertTrue(config["config"]["parameters"]["annotations"])

            login_uri = config["config"]["authentication"]["login"]
            logout_uri = config["config"]["authentication"]["logout"]

            self.assertEqual(login_uri, "/login")
            self.assertEqual(logout_uri, "/logout")

            response = session.get(f"{server}/{login_uri}")
            # check that the login redirect worked
            self.assertEqual(response.history[0].status_code, 302)
            self.assertEqual(response.url, f"{server}/")

            config = session.get(f"{server}/api/v0.2/config").json()
            userinfo = session.get(f"{server}/api/v0.2/userinfo").json()
            self.assertTrue(userinfo["userinfo"]["is_authenticated"])
            self.assertEqual(userinfo["userinfo"]["username"], "test_account")
            self.assertTrue(config["config"]["parameters"]["annotations"])

            response = session.get(f"{server}/{logout_uri}")
            # check that the logout redirect worked
            self.assertEqual(response.history[0].status_code, 302)
            self.assertEqual(response.url, f"{server}/")
            config = session.get(f"{server}/api/v0.2/config").json()
            userinfo = session.get(f"{server}/api/v0.2/userinfo").json()
            self.assertFalse(userinfo["userinfo"]["is_authenticated"])
            self.assertIsNone(userinfo["userinfo"]["username"])
            self.assertTrue(config["config"]["parameters"]["annotations"])
