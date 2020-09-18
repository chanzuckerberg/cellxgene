import os
import unittest
from unittest import mock
from unittest.mock import patch
import tempfile

import requests

from server.common.config.app_config import AppConfig
from server.common.errors import ConfigurationError
from server.common.utils.utils import find_available_port
from server.test import PROJECT_ROOT, test_server, FIXTURES_ROOT


# NOTE, there are more tests that should be written for AppConfig.
# this is just a start.

def mockenv(**envvars):
    return mock.patch.dict(os.environ, envvars)


class AppConfigTest(unittest.TestCase):
    def test_update(self):
        config = AppConfig()
        config.update_server_config(app__verbose=True, multi_dataset__dataroot="datadir")
        vars = config.server_config.changes_from_default()
        self.assertCountEqual(vars, [("app__verbose", True, False), ("multi_dataset__dataroot", "datadir", None)])

        config = AppConfig()
        config.update_default_dataset_config(app__scripts=(), app__inline_scripts=())
        vars = config.server_config.changes_from_default()
        self.assertCountEqual(vars, [])

        config = AppConfig()
        config.update_default_dataset_config(app__scripts=[], app__inline_scripts=[])
        vars = config.default_dataset_config.changes_from_default()
        self.assertCountEqual(vars, [])

        config = AppConfig()
        config.update_default_dataset_config(app__scripts=("a", "b"), app__inline_scripts=["c", "d"])
        vars = config.default_dataset_config.changes_from_default()
        self.assertCountEqual(vars, [("app__scripts", ["a", "b"], []), ("app__inline_scripts", ["c", "d"], [])])

    def test_multi_dataset(self):

        config = AppConfig()
        # test for illegal url_dataroots
        for illegal in ("../b", "!$*", "\\n", "", "(bad)"):
            config.update_server_config(
                multi_dataset__dataroot={"tag": {"base_url": illegal, "dataroot": "{PROJECT_ROOT}/example-dataset"}}
            )
            with self.assertRaises(ConfigurationError):
                config.complete_config()

        # test for legal url_dataroots
        for legal in ("d", "this.is-okay_", "a/b"):
            config.update_server_config(
                multi_dataset__dataroot={"tag": {"base_url": legal, "dataroot": "{PROJECT_ROOT}/example-dataset"}}
            )
            config.complete_config()

        # test that multi dataroots work end to end
        config.update_server_config(
            multi_dataset__dataroot=dict(
                s1=dict(dataroot=f"{PROJECT_ROOT}/example-dataset", base_url="set1/1/2"),
                s2=dict(dataroot=f"{FIXTURES_ROOT}", base_url="set2"),
                s3=dict(dataroot=f"{FIXTURES_ROOT}", base_url="set3"),
            )
        )

        # Change this default to test if the dataroot overrides below work.
        config.update_default_dataset_config(app__about_legal_tos="tos_default.html")

        # specialize the configs for set1
        config.add_dataroot_config(
            "s1", user_annotations__enable=False, diffexp__enable=True, app__about_legal_tos="tos_set1.html"
        )

        # specialize the configs for set2
        config.add_dataroot_config(
            "s2", user_annotations__enable=True, diffexp__enable=False, app__about_legal_tos="tos_set2.html"
        )

        # no specializations for set3 (they get the default dataset config)
        config.complete_config()

        with test_server(app_config=config) as server:
            session = requests.Session()

            response = session.get(f"{server}/set1/1/2/pbmc3k.h5ad/api/v0.2/config")
            data_config = response.json()
            assert data_config["config"]["displayNames"]["dataset"] == "pbmc3k"
            assert data_config["config"]["parameters"]["annotations"] is False
            assert data_config["config"]["parameters"]["disable-diffexp"] is False
            assert data_config["config"]["parameters"]["about_legal_tos"] == "tos_set1.html"

            response = session.get(f"{server}/set2/pbmc3k.cxg/api/v0.2/config")
            data_config = response.json()
            assert data_config["config"]["displayNames"]["dataset"] == "pbmc3k"
            assert data_config["config"]["parameters"]["annotations"] is True
            assert data_config["config"]["parameters"]["about_legal_tos"] == "tos_set2.html"

            response = session.get(f"{server}/set3/pbmc3k.cxg/api/v0.2/config")
            data_config = response.json()
            assert data_config["config"]["displayNames"]["dataset"] == "pbmc3k"
            assert data_config["config"]["parameters"]["annotations"] is True
            assert data_config["config"]["parameters"]["disable-diffexp"] is False
            assert data_config["config"]["parameters"]["about_legal_tos"] == "tos_default.html"

            response = session.get(f"{server}/health")
            assert response.json()["status"] == "pass"

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

    def test_api_base_url(self):

        # test the api_base_url feature, and that it can contain a path
        config = AppConfig()
        backend_port = find_available_port("localhost", 10000)
        config.update_server_config(
            app__api_base_url=f"http://localhost:{backend_port}/additional/path",
            multi_dataset__dataroot=f"{PROJECT_ROOT}/example-dataset"
        )

        config.complete_config()

        with test_server(["-p", str(backend_port)], app_config=config) as server:
            session = requests.Session()
            self.assertEqual(server, f"http://localhost:{backend_port}")
            response = session.get(f"{server}/additional/path/d/pbmc3k.h5ad/api/v0.2/config")
            self.assertEqual(response.status_code, 200)
            data_config = response.json()
            self.assertEqual(data_config["config"]["displayNames"]["dataset"], "pbmc3k")

            # test the health check at the correct url
            response = session.get(f"{server}/additional/path/health")
            assert response.json()["status"] == "pass"

            # also check that the old URL still works.
            # NOTE:  this old URL location will soon be deprecated, and when that happens
            # this check can be removed.
            response = session.get(f"{server}/health")
            assert response.json()["status"] == "pass"

    def test_configfile_with_specialization(self):
        # test that per_dataset_config config load the default config, then the specialized config

        with tempfile.TemporaryDirectory() as tempdir:
            configfile = os.path.join(tempdir, "config.yaml")
            with open(configfile, "w") as fconfig:
                config = """
                server:
                    multi_dataset:
                        dataroot:
                            test:
                                base_url: test
                                dataroot: fake_dataroot

                dataset:
                    user_annotations:
                        enable: false
                        type: hosted_tiledb_array
                        hosted_tiledb_array:
                            db_uri: fake_db_uri
                            hosted_file_directory: fake_dir

                per_dataset_config:
                    test:
                        user_annotations:
                            enable: true
                """
                fconfig.write(config)

            app_config = AppConfig()
            app_config.update_from_config_file(configfile)

            test_config = app_config.dataroot_config["test"]

            # test config from default
            self.assertEqual(test_config.user_annotations__type, "hosted_tiledb_array")
            self.assertEqual(test_config.user_annotations__hosted_tiledb_array__db_uri, "fake_db_uri")

            # test config from specialization
            self.assertTrue(test_config.user_annotations__enable)
