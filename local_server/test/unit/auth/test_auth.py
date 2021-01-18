import unittest

import requests

from local_server.common.config.app_config import AppConfig
from local_server.test import FIXTURES_ROOT, test_server


class AuthTest(unittest.TestCase):
    def setUp(self):
        self.dataset_dataroot = FIXTURES_ROOT

    def test_auth_none(self):
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(authentication__type=None, multi_dataset__dataroot=self.dataset_dataroot)
        app_config.update_default_dataset_config(user_annotations__enable=False)

        app_config.complete_config()

        with test_server(app_config=app_config) as server:
            session = requests.Session()
            config = session.get(f"{server}/d/pbmc3k-CSC-gz.h5ad/api/v0.2/config").json()
            userinfo = session.get(f"{server}/d/pbmc3k-CSC-gz.h5ad/api/v0.2/userinfo").json()
            self.assertNotIn("authentication", config["config"])
            self.assertIsNone(userinfo)

    def test_auth_session(self):
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(authentication__type="session", multi_dataset__dataroot=self.dataset_dataroot)
        app_config.update_default_dataset_config(user_annotations__enable=True)
        app_config.complete_config()

        with test_server(app_config=app_config) as server:
            session = requests.Session()
            config = session.get(f"{server}/d/pbmc3k-CSC-gz.h5ad/api/v0.2/config").json()
            userinfo = session.get(f"{server}/d/pbmc3k-CSC-gz.h5ad/api/v0.2/userinfo").json()

            self.assertFalse(config["config"]["authentication"]["requires_client_login"])
            self.assertTrue(userinfo["userinfo"]["is_authenticated"])
            self.assertEqual(userinfo["userinfo"]["username"], "anonymous")

    def test_auth_test(self):
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(authentication__type="test")
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

        with test_server(app_config=app_config) as server:
            session = requests.Session()

            # auth datasets
            config = session.get(f"{server}/auth/pbmc3k-CSC-gz.h5ad/api/v0.2/config").json()
            userinfo = session.get(f"{server}/auth/pbmc3k-CSC-gz.h5ad/api/v0.2/userinfo").json()

            self.assertFalse(userinfo["userinfo"]["is_authenticated"])
            self.assertIsNone(userinfo["userinfo"]["username"])
            self.assertTrue(config["config"]["authentication"]["requires_client_login"])
            self.assertTrue(config["config"]["parameters"]["annotations"])

            login_uri = config["config"]["authentication"]["login"]
            logout_uri = config["config"]["authentication"]["logout"]

            self.assertEqual(login_uri, "/login?dataset=auth/pbmc3k-CSC-gz.h5ad")
            self.assertEqual(logout_uri, "/logout?dataset=auth/pbmc3k-CSC-gz.h5ad")

            r = session.get(f"{server}/{login_uri}")
            # check that the login redirect worked
            self.assertEqual(r.history[0].status_code, 302)
            self.assertEqual(r.url, f"{server}/auth/pbmc3k-CSC-gz.h5ad")

            config = session.get(f"{server}/auth/pbmc3k-CSC-gz.h5ad/api/v0.2/config").json()
            userinfo = session.get(f"{server}/auth/pbmc3k-CSC-gz.h5ad/api/v0.2/userinfo").json()
            self.assertTrue(userinfo["userinfo"]["is_authenticated"])
            self.assertEqual(userinfo["userinfo"]["username"], "test_account")
            self.assertEqual(userinfo["userinfo"]["picture"], None)
            self.assertTrue(config["config"]["parameters"]["annotations"])

            r = session.get(f"{server}/{logout_uri}")
            # check that the logout redirect worked
            self.assertEqual(r.history[0].status_code, 302)
            self.assertEqual(r.url, f"{server}/auth/pbmc3k-CSC-gz.h5ad")
            config = session.get(f"{server}/auth/pbmc3k-CSC-gz.h5ad/api/v0.2/config").json()
            userinfo = session.get(f"{server}/auth/pbmc3k-CSC-gz.h5ad/api/v0.2/userinfo").json()
            self.assertFalse(userinfo["userinfo"]["is_authenticated"])
            self.assertIsNone(userinfo["userinfo"]["username"])
            self.assertTrue(config["config"]["parameters"]["annotations"])

            # no-auth datasets
            config = session.get(f"{server}/no-auth/pbmc3k-CSC-gz.h5ad/api/v0.2/config").json()
            userinfo = session.get(f"{server}/no-auth/pbmc3k-CSC-gz.h5ad/api/v0.2/userinfo").json()
            self.assertIsNone(userinfo)
            self.assertFalse(config["config"]["parameters"]["annotations"])

            # login with a picture
            session.get(f"{server}/{login_uri}&picture=myimage.png")
            userinfo = session.get(f"{server}/auth/pbmc3k-CSC-gz.h5ad/api/v0.2/userinfo").json()
            self.assertTrue(userinfo["userinfo"]["is_authenticated"])
            self.assertEqual(userinfo["userinfo"]["picture"], "myimage.png")

    def test_auth_test_single(self):
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(
            authentication__type="test", single_dataset__datapath=f"{self.dataset_dataroot}/pbmc3k-CSC-gz.h5ad"
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
