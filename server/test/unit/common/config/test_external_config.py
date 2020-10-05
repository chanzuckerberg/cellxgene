import os
import tempfile
import unittest
from unittest.mock import patch

import requests

from server.common.config.app_config import ExternalConfig
from server.common.config.app_config import AppConfig
from server.test import test_server, FIXTURES_ROOT


class TestExternalConfig(unittest.TestCase):
    def test_type_convert(self):
        # The values from environment variables and aws secrets are returned as strings.
        # These values need to be converted to the proper types.

        self.assertEqual(ExternalConfig.convert_value("1"), int(1))
        self.assertEqual(ExternalConfig.convert_value("1.1"), float(1.1))
        self.assertEqual(ExternalConfig.convert_value("string"), "string")
        self.assertEqual(ExternalConfig.convert_value("true"), True)
        self.assertEqual(ExternalConfig.convert_value("True"), True)
        self.assertEqual(ExternalConfig.convert_value("false"), False)
        self.assertEqual(ExternalConfig.convert_value("False"), False)
        self.assertEqual(ExternalConfig.convert_value("False"), False)
        self.assertEqual(ExternalConfig.convert_value("null"), None)
        self.assertEqual(ExternalConfig.convert_value("None"), None)
        self.assertEqual(ExternalConfig.convert_value("{'a':10, 'b':'string'}"), dict(a=int(10), b="string"))

    def test_environment_variable(self):
        with tempfile.TemporaryDirectory() as tempdir:
            configfile = os.path.join(tempdir, "config.yaml")
            with open(configfile, "w") as fconfig:
                config = """
                external:
                  environment:
                    - name: DATAPATH
                      path: [server, single_dataset, datapath]
                      required: true
                    - name: DIFFEXP
                      path: [dataset, diffexp, enable]
                      required: true
                """
                fconfig.write(config)

            env = os.environ
            env["DATAPATH"] = f"{FIXTURES_ROOT}/pbmc3k.cxg"
            env["DIFFEXP"] = "False"
            with test_server(command_line_args=["-c", configfile], env=env) as server:
                session = requests.Session()
                response = session.get(f"{server}/api/v0.2/config")
                data_config = response.json()
                self.assertEqual(data_config["config"]["displayNames"]["dataset"], "pbmc3k")
                self.assertTrue(data_config["config"]["parameters"]["disable-diffexp"])

            env["DATAPATH"] = f"{FIXTURES_ROOT}/a95c59b4-7f5d-4b80-ad53-a694834ca18b.h5ad"
            env["DIFFEXP"] = "True"
            with test_server(command_line_args=["-c", configfile], env=env) as server:
                session = requests.Session()
                response = session.get(f"{server}/api/v0.2/config")
                data_config = response.json()
                self.assertEqual(
                    data_config["config"]["displayNames"]["dataset"], "a95c59b4-7f5d-4b80-ad53-a694834ca18b"
                )
                self.assertFalse(data_config["config"]["parameters"]["disable-diffexp"])

    @patch("server.common.config.external_config.get_secret_key")
    def test_aws_secrets_manager(self, mock_get_secret_key):
        mock_get_secret_key.return_value = {
            "oauth_client_secret": "mock_oauth_secret",
            "db_uri": "mock_db_uri",
        }
        with tempfile.TemporaryDirectory() as tempdir:
            configfile = os.path.join(tempdir, "config.yaml")
            with open(configfile, "w") as fconfig:
                config = """
                external:
                  aws_secrets_manager:
                    region: us-west-2
                    secrets:
                      - name: my_secret
                        values:
                          - key: flask_secret_key
                            path: [server, app, flask_secret_key]
                            required: false
                          - key: db_uri
                            path: [dataset, user_annotations, hosted_tiledb_array, db_uri]
                            required: true
                          - key: oauth_client_secret
                            path: [server, authentication, params_oauth, client_secret]
                            required: true
                """
                fconfig.write(config)

            app_config = AppConfig()
            app_config.update_from_config_file(configfile)
            app_config.server_config.single_dataset__datapath = f"{FIXTURES_ROOT}/pbmc3k.cxg"
            app_config.server_config.app__flask_secret_key = "original"
            app_config.server_config.single_dataset__datapath = f"{FIXTURES_ROOT}/pbmc3k.cxg"

            app_config.complete_config()

            self.assertEqual(app_config.server_config.app__flask_secret_key, "original")
            self.assertEqual(app_config.server_config.authentication__params_oauth__client_secret, "mock_oauth_secret")
            self.assertEqual(
                app_config.default_dataset_config.user_annotations__hosted_tiledb_array__db_uri, "mock_db_uri"
            )
