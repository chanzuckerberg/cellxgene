import unittest
from server.common.app_config import AppConfig
from server.test import PROJECT_ROOT, test_server
import requests


class AuthTest(unittest.TestCase):
    def test_auth_none(self):
        c = AppConfig()
        c.update_server_config(
            authentication__type=None, multi_dataset__dataroot=f"{PROJECT_ROOT}/server/test/test_datasets"
        )
        c.update_default_dataset_config(user_annotations__enable=False)

        c.complete_config()

        with test_server(app_config=c) as server:
            session = requests.Session()
            r = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/config")
            data_config = r.json()
            assert "authentication" not in data_config["config"]

    def test_auth_session(self):
        c = AppConfig()
        c.update_server_config(
            authentication__type="session", multi_dataset__dataroot=f"{PROJECT_ROOT}/server/test/test_datasets"
        )
        c.update_default_dataset_config(user_annotations__enable=True)
        c.complete_config()

        with test_server(app_config=c) as server:
            session = requests.Session()
            r = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/config")
            data_config = r.json()
            assert data_config["config"]["authentication"]["is_authenticated"]
            assert not data_config["config"]["authentication"]["requires_client_login"]
            assert data_config["config"]["authentication"]["username"] == "anonymous"

    def test_auth_test(self):
        c = AppConfig()
        c.update_server_config(authentication__type="test")
        c.update_server_config(
            multi_dataset__dataroot=dict(
                a1=dict(dataroot=f"{PROJECT_ROOT}/server/test/test_datasets", base_url="auth"),
                a2=dict(dataroot=f"{PROJECT_ROOT}/server/test/test_datasets", base_url="no-auth"),
            )
        )

        # specialize the configs
        c.add_dataroot_config("a1", app__authentication_enable=True, user_annotations__enable=True)
        c.add_dataroot_config("a2", app__authentication_enable=False, user_annotations__enable=False)

        c.complete_config()

        with test_server(app_config=c) as server:
            session = requests.Session()

            # auth datasets
            r = session.get(f"{server}/auth/pbmc3k.cxg/api/v0.2/config")
            data_config = r.json()
            assert not data_config["config"]["authentication"]["is_authenticated"]
            assert data_config["config"]["authentication"]["requires_client_login"]
            assert data_config["config"]["authentication"]["username"] is None
            assert data_config["config"]["parameters"]["annotations"]

            login_uri = data_config["config"]["authentication"]["login"]
            logout_uri = data_config["config"]["authentication"]["logout"]

            assert login_uri == "/login?dataset=auth/pbmc3k.cxg"
            assert logout_uri == "/logout?dataset=auth/pbmc3k.cxg"

            r = session.get(f"{server}/{login_uri}")
            # check that the login redirect worked
            assert r.history[0].status_code == 302
            assert r.url == f"{server}/auth/pbmc3k.cxg/"

            r = session.get(f"{server}/auth/pbmc3k.cxg/api/v0.2/config")
            data_config = r.json()
            assert data_config["config"]["authentication"]["is_authenticated"]
            assert data_config["config"]["authentication"]["username"] == "test_account"
            assert data_config["config"]["parameters"]["annotations"]

            r = session.get(f"{server}/{logout_uri}")
            # check that the logout redirect worked
            assert r.history[0].status_code == 302
            assert r.url == f"{server}/auth/pbmc3k.cxg/"
            r = session.get(f"{server}/auth/pbmc3k.cxg/api/v0.2/config")
            data_config = r.json()
            assert not data_config["config"]["authentication"]["is_authenticated"]
            assert data_config["config"]["authentication"]["username"] is None
            assert data_config["config"]["parameters"]["annotations"]

            # no-auth datasets
            r = session.get(f"{server}/no-auth/pbmc3k.cxg/api/v0.2/config")
            data_config = r.json()
            assert "authentication" not in data_config["config"]
            assert not data_config["config"]["parameters"]["annotations"]

    def test_auth_test_single(self):
        c = AppConfig()
        c.update_server_config(
            authentication__type="test",
            single_dataset__datapath=f"{PROJECT_ROOT}/server/test/test_datasets/pbmc3k.cxg")

        c.complete_config()

        with test_server(app_config=c) as server:
            session = requests.Session()

            r = session.get(f"{server}/api/v0.2/config")
            data_config = r.json()
            assert not data_config["config"]["authentication"]["is_authenticated"]
            assert data_config["config"]["authentication"]["requires_client_login"]
            assert data_config["config"]["authentication"]["username"] is None
            assert data_config["config"]["parameters"]["annotations"]

            login_uri = data_config["config"]["authentication"]["login"]
            logout_uri = data_config["config"]["authentication"]["logout"]

            assert login_uri == "/login"
            assert logout_uri == "/logout"

            r = session.get(f"{server}/{login_uri}")
            # check that the login redirect worked
            assert r.history[0].status_code == 302
            assert r.url == f"{server}/"

            r = session.get(f"{server}/api/v0.2/config")
            data_config = r.json()
            assert data_config["config"]["authentication"]["is_authenticated"]
            assert data_config["config"]["authentication"]["username"] == "test_account"
            assert data_config["config"]["parameters"]["annotations"]

            r = session.get(f"{server}/{logout_uri}")
            # check that the logout redirect worked
            assert r.history[0].status_code == 302
            assert r.url == f"{server}/"
            r = session.get(f"{server}/api/v0.2/config")
            data_config = r.json()
            assert not data_config["config"]["authentication"]["is_authenticated"]
            assert data_config["config"]["authentication"]["username"] is None
            assert data_config["config"]["parameters"]["annotations"]
