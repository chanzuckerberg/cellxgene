import os
import unittest
from unittest import mock
from unittest.mock import patch

import requests

from server.common.app_config import AppConfig
from server.common.errors import ConfigurationError
from server.common.utils.utils import find_available_port
from server.test import PROJECT_ROOT, test_server, FIXTURES_ROOT


# NOTE, there are more tests that should be written for AppConfig.
# this is just a start.

def mockenv(**envvars):
    return mock.patch.dict(os.environ, envvars)


class AppConfigTest(unittest.TestCase):
    def test_update(self):
        c = AppConfig()
        c.update_server_config(app__verbose=True, multi_dataset__dataroot="datadir")
        v = c.server_config.changes_from_default()
        self.assertCountEqual(v, [("app__verbose", True, False), ("multi_dataset__dataroot", "datadir", None)])

        c = AppConfig()
        c.update_default_dataset_config(app__scripts=(), app__inline_scripts=())
        v = c.server_config.changes_from_default()
        self.assertCountEqual(v, [])

        c = AppConfig()
        c.update_default_dataset_config(app__scripts=[], app__inline_scripts=[])
        v = c.default_dataset_config.changes_from_default()
        self.assertCountEqual(v, [])

        c = AppConfig()
        c.update_default_dataset_config(app__scripts=("a", "b"), app__inline_scripts=["c", "d"])
        v = c.default_dataset_config.changes_from_default()
        self.assertCountEqual(v, [("app__scripts", ["a", "b"], []), ("app__inline_scripts", ["c", "d"], [])])

    def test_multi_dataset(self):

        c = AppConfig()
        # test for illegal url_dataroots
        for illegal in ("../b", "!$*", "\\n", "", "(bad)"):
            c.update_server_config(
                multi_dataset__dataroot={"tag": {"base_url": illegal, "dataroot": "{PROJECT_ROOT}/example-dataset"}}
            )
            with self.assertRaises(ConfigurationError):
                c.complete_config()

        # test for legal url_dataroots
        for legal in ("d", "this.is-okay_", "a/b"):
            c.update_server_config(
                multi_dataset__dataroot={"tag": {"base_url": legal, "dataroot": "{PROJECT_ROOT}/example-dataset"}}
            )
            c.complete_config()

        # test that multi dataroots work end to end
        c.update_server_config(
            multi_dataset__dataroot=dict(
                s1=dict(dataroot=f"{PROJECT_ROOT}/example-dataset", base_url="set1/1/2"),
                s2=dict(dataroot=f"{FIXTURES_ROOT}", base_url="set2"),
                s3=dict(dataroot=f"{FIXTURES_ROOT}", base_url="set3"),
            )
        )

        # Change this default to test if the dataroot overrides below work.
        c.update_default_dataset_config(app__about_legal_tos="tos_default.html")

        # specialize the configs for set1
        c.add_dataroot_config(
            "s1", user_annotations__enable=False, diffexp__enable=True, app__about_legal_tos="tos_set1.html"
        )

        # specialize the configs for set2
        c.add_dataroot_config(
            "s2", user_annotations__enable=True, diffexp__enable=False, app__about_legal_tos="tos_set2.html"
        )

        # no specializations for set3 (they get the default dataset config)
        c.complete_config()

        with test_server(app_config=c) as server:
            session = requests.Session()

            r = session.get(f"{server}/set1/1/2/pbmc3k.h5ad/api/v0.2/config")
            data_config = r.json()
            assert data_config["config"]["displayNames"]["dataset"] == "pbmc3k"
            assert data_config["config"]["parameters"]["annotations"] is False
            assert data_config["config"]["parameters"]["disable-diffexp"] is False
            assert data_config["config"]["parameters"]["about_legal_tos"] == "tos_set1.html"

            r = session.get(f"{server}/set2/pbmc3k.cxg/api/v0.2/config")
            data_config = r.json()
            assert data_config["config"]["displayNames"]["dataset"] == "pbmc3k"
            assert data_config["config"]["parameters"]["annotations"] is True
            assert data_config["config"]["parameters"]["about_legal_tos"] == "tos_set2.html"

            r = session.get(f"{server}/set3/pbmc3k.cxg/api/v0.2/config")
            data_config = r.json()
            assert data_config["config"]["displayNames"]["dataset"] == "pbmc3k"
            assert data_config["config"]["parameters"]["annotations"] is True
            assert data_config["config"]["parameters"]["disable-diffexp"] is False
            assert data_config["config"]["parameters"]["about_legal_tos"] == "tos_default.html"

            r = session.get(f"{server}/health")
            assert r.json()["status"] == "pass"

    @mockenv(CXG_AWS_SECRET_NAME="TESTING", CXG_AWS_SECRET_REGION_NAME="TEST_REGION")
    @patch('server.common.aws_secret_utils.get_secret_key')
    def test_get_config_vars_from_aws_secrets(self, mock_get_secret_key):
        mock_get_secret_key.return_value = {
            "flask_secret_key": "mock_flask_secret",
            "oauth_client_secret": "mock_oauth_secret",
            "db_uri": "mock_db_uri"
        }

        config = AppConfig()

        with self.assertLogs(level="INFO") as logger:
            from server.common.aws_secret_utils import handle_config_from_secret
            # should not throw error
            # "AttributeError: 'XConfig' object has no attribute 'x'"
            handle_config_from_secret(config)

            # should log 3 lines (one for each var set from a secret)
            self.assertEqual(len(logger.output), 3)
            self.assertIn('INFO:root:set app__flask_secret_key from secret', logger.output[0])
            self.assertIn('INFO:root:set authentication__params_oauth__client_secret from secret', logger.output[1])
            self.assertIn('INFO:root:set user_annotations__hosted_tiledb_array__db_uri from secret', logger.output[2])
            self.assertEqual(config.server_config.app__flask_secret_key, "mock_flask_secret")
            self.assertEqual(config.server_config.authentication__params_oauth__client_secret, "mock_oauth_secret")
            self.assertEqual(config.default_dataset_config.user_annotations__hosted_tiledb_array__db_uri, "mock_db_uri")

    def test_backend_base_url(self):

        # test the backend_base_url feature, and that it can contain a path
        c = AppConfig()
        backend_port = find_available_port("localhost", 10000)
        c.update_server_config(
            app__backend_base_url=f"http://localhost:{backend_port}/additional/path/before/dataroot",
            multi_dataset__dataroot=f"{PROJECT_ROOT}/example-dataset"
        )

        c.complete_config()

        with test_server(["-p", str(backend_port)], app_config=c) as server:
            session = requests.Session()
            self.assertEqual(server, f"http://localhost:{backend_port}")
            response = session.get(f"{server}/additional/path/before/dataroot/d/pbmc3k.h5ad/api/v0.2/config")
            self.assertEqual(response.status_code, 200)
            data_config = response.json()
            self.assertEqual(data_config["config"]["displayNames"]["dataset"], "pbmc3k")
