import os
from unittest.mock import patch

import requests

from server.common.errors import ConfigurationError
from server.common.config.app_config import AppConfig
from test.unit import test_server
from test import FIXTURES_ROOT
from server.common.utils.type_conversion_utils import convert_string_to_value
from test.unit.common.config import ConfigTests


class TestExternalConfig(ConfigTests):
    def test_type_convert(self):
        # The values from environment variables are returned as strings.
        # These values need to be converted to the proper types.

        self.assertEqual(convert_string_to_value("1"), int(1))
        self.assertEqual(convert_string_to_value("1.1"), float(1.1))
        self.assertEqual(convert_string_to_value("string"), "string")
        self.assertEqual(convert_string_to_value("true"), True)
        self.assertEqual(convert_string_to_value("True"), True)
        self.assertEqual(convert_string_to_value("false"), False)
        self.assertEqual(convert_string_to_value("False"), False)
        self.assertEqual(convert_string_to_value("null"), None)
        self.assertEqual(convert_string_to_value("None"), None)
        self.assertEqual(convert_string_to_value("{'a':10, 'b':'string'}"), dict(a=int(10), b="string"))

    def test_environment_variable(self):
        configfile = self.custom_external_config(
            environment=[
                dict(name="DATAPATH", path=["server", "single_dataset", "datapath"], required=True),
                dict(name="DIFFEXP", path=["dataset", "diffexp", "enable"], required=True),
            ],
            config_file_name="environment_external_config.yaml",
        )

        env = os.environ
        env["DATAPATH"] = f"{FIXTURES_ROOT}/pbmc3k-CSC-gz.h5ad"
        env["DIFFEXP"] = "False"
        with test_server(command_line_args=["-c", configfile], env=env) as server:
            session = requests.Session()
            response = session.get(f"{server}/api/v0.2/config")
            data_config = response.json()
            self.assertEqual(data_config["config"]["displayNames"]["dataset"], "pbmc3k-CSC-gz")
            self.assertTrue(data_config["config"]["parameters"]["disable-diffexp"])

        env["DATAPATH"] = f"{FIXTURES_ROOT}/a95c59b4-7f5d-4b80-ad53-a694834ca18b.h5ad"
        env["DIFFEXP"] = "True"
        with test_server(command_line_args=["-c", configfile], env=env) as server:
            session = requests.Session()
            response = session.get(f"{server}/api/v0.2/config")
            data_config = response.json()
            self.assertEqual(data_config["config"]["displayNames"]["dataset"], "a95c59b4-7f5d-4b80-ad53-a694834ca18b")
            self.assertFalse(data_config["config"]["parameters"]["disable-diffexp"])

    def test_environment_variable_errors(self):

        # no name
        app_config = AppConfig()
        app_config.external_config.environment = [dict(required=True, path=["this", "is", "a", "path"])]
        with self.assertRaises(ConfigurationError) as config_error:
            app_config.complete_config()
        self.assertEqual(config_error.exception.message, "environment: 'name' is missing")

        # required has wrong type
        app_config = AppConfig()
        app_config.external_config.environment = [
            dict(name="myenvar", required="optional", path=["this", "is", "a", "path"])
        ]
        with self.assertRaises(ConfigurationError) as config_error:
            app_config.complete_config()
        self.assertEqual(config_error.exception.message, "environment: 'required' must be a bool")

        # no path
        app_config = AppConfig()
        app_config.external_config.environment = [dict(name="myenvar", required=True)]
        with self.assertRaises(ConfigurationError) as config_error:
            app_config.complete_config()
        self.assertEqual(config_error.exception.message, "environment: 'path' is missing")

        # required environment variable is not set
        app_config = AppConfig()
        app_config.external_config.environment = [
            dict(name="THIS_ENV_IS_NOT_SET", required=True, path=["this", "is", "a", "path"])
        ]
        with self.assertRaises(ConfigurationError) as config_error:
            app_config.complete_config()
        self.assertEqual(config_error.exception.message, "required environment variable 'THIS_ENV_IS_NOT_SET' not set")
<<<<<<< HEAD:backend/test/test_server/unit/common/config/test_external_config.py
=======

    @patch("server.common.config.external_config.get_secret_key")
    def test_aws_secrets_manager(self, mock_get_secret_key):
        mock_get_secret_key.return_value = {
            "flask_secret_key": "mock_flask_secret_key",
        }
        configfile = self.custom_external_config(
            aws_secrets_manager_region="us-west-2",
            aws_secrets_manager_secrets=[
                dict(
                    name="my_secret",
                    values=[
                        dict(key="flask_secret_key", path=["server", "app", "flask_secret_key"], required=True),
                    ],
                )
            ],
            config_file_name="secret_external_config.yaml",
        )

        app_config = AppConfig()
        app_config.update_from_config_file(configfile)
        app_config.server_config.single_dataset__datapath = f"{FIXTURES_ROOT}/pbmc3k-CSC-gz.h5ad"

        app_config.complete_config()

        self.assertEqual(app_config.server_config.app__flask_secret_key, "mock_flask_secret_key")

    @patch("server.common.config.external_config.get_secret_key")
    def test_aws_secrets_manager_error(self, mock_get_secret_key):
        mock_get_secret_key.return_value = {
            "db_uri": "mock_db_uri",
        }

        # no region
        app_config = AppConfig()
        app_config.external_config.aws_secrets_manager__region = None
        app_config.external_config.aws_secrets_manager__secrets = [
            dict(name="secret1", values=[dict(key="key1", required=True, path=["this", "is", "my", "path"])])
        ]
        with self.assertRaises(ConfigurationError) as config_error:
            app_config.complete_config()
        self.assertEqual(
            config_error.exception.message,
            "Invalid type for attribute: aws_secrets_manager__region, expected type str, got NoneType",
        )

        # missing secret name
        app_config = AppConfig()
        app_config.external_config.aws_secrets_manager__region = "us-west-2"
        app_config.external_config.aws_secrets_manager__secrets = [
            dict(values=[dict(key="db_uri", required=True, path=["this", "is", "my", "path"])])
        ]
        with self.assertRaises(ConfigurationError) as config_error:
            app_config.complete_config()
        self.assertEqual(config_error.exception.message, "aws_secrets_manager: 'name' is missing")

        # secret name wrong type
        app_config = AppConfig()
        app_config.external_config.aws_secrets_manager__region = "us-west-2"
        app_config.external_config.aws_secrets_manager__secrets = [
            dict(name=1, values=[dict(key="db_uri", required=True, path=["this", "is", "my", "path"])])
        ]
        with self.assertRaises(ConfigurationError) as config_error:
            app_config.complete_config()
        self.assertEqual(config_error.exception.message, "aws_secrets_manager: 'name' must be a string")

        # missing values name
        app_config = AppConfig()
        app_config.external_config.aws_secrets_manager__region = "us-west-2"
        app_config.external_config.aws_secrets_manager__secrets = [dict(name="mysecret")]
        with self.assertRaises(ConfigurationError) as config_error:
            app_config.complete_config()
        self.assertEqual(config_error.exception.message, "aws_secrets_manager: 'values' is missing")

        # values wrong type
        app_config = AppConfig()
        app_config.external_config.aws_secrets_manager__region = "us-west-2"
        app_config.external_config.aws_secrets_manager__secrets = [
            dict(name="mysecret", values=dict(key="db_uri", required=True, path=["this", "is", "my", "path"]))
        ]
        with self.assertRaises(ConfigurationError) as config_error:
            app_config.complete_config()
        self.assertEqual(config_error.exception.message, "aws_secrets_manager: 'values' must be a list")

        # entry missing key
        app_config = AppConfig()
        app_config.external_config.aws_secrets_manager__region = "us-west-2"
        app_config.external_config.aws_secrets_manager__secrets = [
            dict(name="mysecret", values=[dict(required=True, path=["this", "is", "my", "path"])])
        ]
        with self.assertRaises(ConfigurationError) as config_error:
            app_config.complete_config()
        self.assertEqual(config_error.exception.message, "missing 'key' in secret values: mysecret")

        # entry required is wrong type
        app_config = AppConfig()
        app_config.external_config.aws_secrets_manager__region = "us-west-2"
        app_config.external_config.aws_secrets_manager__secrets = [
            dict(name="mysecret", values=[dict(key="db_uri", required="optional", path=["this", "is", "my", "path"])])
        ]
        with self.assertRaises(ConfigurationError) as config_error:
            app_config.complete_config()
        self.assertEqual(config_error.exception.message, "wrong type for 'required' in secret values: mysecret")

        # entry missing path
        app_config = AppConfig()
        app_config.external_config.aws_secrets_manager__region = "us-west-2"
        app_config.external_config.aws_secrets_manager__secrets = [
            dict(name="mysecret", values=[dict(key="db_uri", required=True)])
        ]
        with self.assertRaises(ConfigurationError) as config_error:
            app_config.complete_config()
        self.assertEqual(config_error.exception.message, "missing 'path' in secret values: mysecret")

        # secret missing required key
        app_config = AppConfig()
        app_config.external_config.aws_secrets_manager__region = "us-west-2"
        app_config.external_config.aws_secrets_manager__secrets = [
            dict(
                name="mysecret",
                values=[dict(key="KEY_DOES_NOT_EXIST", required=True, path=["this", "is", "a", "path"])],
            )
        ]
        with self.assertRaises(ConfigurationError) as config_error:
            app_config.complete_config()
        self.assertEqual(config_error.exception.message, "required secret 'mysecret:KEY_DOES_NOT_EXIST' not set")
