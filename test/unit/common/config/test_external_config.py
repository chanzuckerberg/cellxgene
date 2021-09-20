import os

import requests

from common.errors import ConfigurationError
from server.common.config.app_config import AppConfig
from common.utils.type_conversion_utils import convert_string_to_value

from test import FIXTURES_ROOT
from test.unit import test_server
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
