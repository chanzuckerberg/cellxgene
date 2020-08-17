import unittest

import requests

from server.common.app_config import AppConfig
from server.test import FIXTURES_ROOT, test_server


class AuthTest(unittest.TestCase):
    def setUp(self):
        self.dataset_dataroot = FIXTURES_ROOT

    def test_auth_none(self):
        c = AppConfig()
        c.update_server_config(
            authentication__type=None, multi_dataset__dataroot=self.dataset_dataroot
        )
        c.update_default_dataset_config(user_annotations__enable=False)

        c.complete_config()

        with test_server(app_config=c) as server:
            session = requests.Session()
            config = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/config").json()
            userinfo = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/userinfo").json()
            self.assertNotIn("authentication", config["config"])
            self.assertIsNone(userinfo)

    def test_auth_session(self):
        c = AppConfig()
        c.update_server_config(
            authentication__type="session", multi_dataset__dataroot=self.dataset_dataroot
        )
        c.update_default_dataset_config(user_annotations__enable=True)
        c.complete_config()

        with test_server(app_config=c) as server:
            session = requests.Session()
            config = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/config").json()
            userinfo = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/userinfo").json()

            self.assertFalse(config["config"]["authentication"]["requires_client_login"])
            self.assertTrue(userinfo["userinfo"]["is_authenticated"])
            self.assertEqual(userinfo["userinfo"]["username"], "anonymous")

    def test_auth_test(self):
        c = AppConfig()
        c.update_server_config(authentication__type="test")
        c.update_server_config(
            multi_dataset__dataroot=dict(
                a1=dict(dataroot=self.dataset_dataroot, base_url="auth"),
                a2=dict(dataroot=self.dataset_dataroot, base_url="no-auth"),
            )
        )

        # specialize the configs
        c.add_dataroot_config("a1", app__authentication_enable=True, user_annotations__enable=True)
        c.add_dataroot_config("a2", app__authentication_enable=False, user_annotations__enable=False)

        c.complete_config()

        with test_server(app_config=c) as server:
            session = requests.Session()

            # auth datasets
            config = session.get(f"{server}/auth/pbmc3k.cxg/api/v0.2/config").json()
            userinfo = session.get(f"{server}/auth/pbmc3k.cxg/api/v0.2/userinfo").json()

            self.assertFalse(userinfo["userinfo"]["is_authenticated"])
            self.assertIsNone(userinfo["userinfo"]["username"])
            self.assertTrue(config["config"]["authentication"]["requires_client_login"])
            self.assertTrue(config["config"]["parameters"]["annotations"])

            login_uri = config["config"]["authentication"]["login"]
            logout_uri = config["config"]["authentication"]["logout"]

            self.assertEqual(login_uri, "/login?dataset=auth/pbmc3k.cxg")
            self.assertEqual(logout_uri, "/logout?dataset=auth/pbmc3k.cxg")

            r = session.get(f"{server}/{login_uri}")
            # check that the login redirect worked
            self.assertEqual(r.history[0].status_code, 302)
            self.assertEqual(r.url, f"{server}/auth/pbmc3k.cxg/")

            config = session.get(f"{server}/auth/pbmc3k.cxg/api/v0.2/config").json()
            userinfo = session.get(f"{server}/auth/pbmc3k.cxg/api/v0.2/userinfo").json()
            self.assertTrue(userinfo["userinfo"]["is_authenticated"])
            self.assertEqual(userinfo["userinfo"]["username"], "test_account")
            self.assertTrue(config["config"]["parameters"]["annotations"])

            r = session.get(f"{server}/{logout_uri}")
            # check that the logout redirect worked
            self.assertEqual(r.history[0].status_code, 302)
            self.assertEqual(r.url, f"{server}/auth/pbmc3k.cxg/")
            config = session.get(f"{server}/auth/pbmc3k.cxg/api/v0.2/config").json()
            userinfo = session.get(f"{server}/auth/pbmc3k.cxg/api/v0.2/userinfo").json()
            self.assertFalse(userinfo["userinfo"]["is_authenticated"])
            self.assertIsNone(userinfo["userinfo"]["username"])
            self.assertTrue(config["config"]["parameters"]["annotations"])

            # no-auth datasets
            config = session.get(f"{server}/no-auth/pbmc3k.cxg/api/v0.2/config").json()
            userinfo = session.get(f"{server}/no-auth/pbmc3k.cxg/api/v0.2/userinfo").json()
            self.assertIsNone(userinfo)
            self.assertFalse(config["config"]["parameters"]["annotations"])

    def test_auth_test_single(self):
        c = AppConfig()
        c.update_server_config(
            authentication__type="test",
            single_dataset__datapath=f"{self.dataset_dataroot}/pbmc3k.cxg")

        c.complete_config()

        with test_server(app_config=c) as server:
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

            r = session.get(f"{server}/{login_uri}")
            # check that the login redirect worked
            self.assertEqual(r.history[0].status_code, 302)
            self.assertEqual(r.url, f"{server}/")

            config = session.get(f"{server}/api/v0.2/config").json()
            userinfo = session.get(f"{server}/api/v0.2/userinfo").json()
            self.assertTrue(userinfo["userinfo"]["is_authenticated"])
            self.assertEqual(userinfo["userinfo"]["username"], "test_account")
            self.assertTrue(config["config"]["parameters"]["annotations"])

            r = session.get(f"{server}/{logout_uri}")
            # check that the logout redirect worked
            self.assertEqual(r.history[0].status_code, 302)
            self.assertEqual(r.url, f"{server}/")
            config = session.get(f"{server}/api/v0.2/config").json()
            userinfo = session.get(f"{server}/api/v0.2/userinfo").json()
            self.assertFalse(userinfo["userinfo"]["is_authenticated"])
            self.assertIsNone(userinfo["userinfo"]["username"])
            self.assertTrue(config["config"]["parameters"]["annotations"])
