import os
import unittest
from unittest import mock
from unittest.mock import patch

from server.common.config.base_config import BaseConfig
from test import H5AD_FIXTURE

from server.common.config.app_config import AppConfig
from server.common.errors import ConfigurationError
from test.unit.common.config import ConfigTests


def mockenv(**envvars):
    return mock.patch.dict(os.environ, envvars)


class TestServerConfig(ConfigTests):
    def setUp(self):
        self.config_file_name = f"{unittest.TestCase.id(self).split('.')[-1]}.yml"
        self.config = AppConfig()
        self.config.update_server_config(app__flask_secret_key="secret")
        self.config.update_server_config(single_dataset__datapath=H5AD_FIXTURE)
        self.server_config = self.config.server_config
        self.config.complete_config()

        message_list = []

        def noop(message):
            message_list.append(message)

        messagefn = noop
        self.context = dict(messagefn=messagefn, messages=message_list)

    def get_config(self, **kwargs):
        file_name = self.custom_app_config(
            dataset_datapath=f"{H5AD_FIXTURE}", config_file_name=self.config_file_name, **kwargs
        )
        config = AppConfig()
        config.update_from_config_file(file_name)
        return config

    def test_init_raises_error_if_default_config_is_invalid(self):
        invalid_config = self.get_config(port="not_valid")
        with self.assertRaises(ConfigurationError):
            invalid_config.complete_config()

    @patch("server.common.config.server_config.BaseConfig.validate_correct_type_of_configuration_attribute")
    def test_complete_config_checks_all_attr(self, mock_check_attrs):
        mock_check_attrs.side_effect = BaseConfig.validate_correct_type_of_configuration_attribute()
        self.server_config.complete_config(self.context)
        self.assertEqual(mock_check_attrs.call_count, 19)

    def test_handle_app__throws_error_if_port_doesnt_exist(self):
        config = self.get_config(port=99999999)
        with self.assertRaises(ConfigurationError):
            config.server_config.handle_app(self.context)

    @patch("server.common.config.server_config.discover_s3_region_name")
    def test_handle_data_locator_works_for_default_types(self, mock_discover_region_name):
        mock_discover_region_name.return_value = None
        # Default config
        self.assertEqual(self.config.server_config.data_locator__s3__region_name, None)
        # hard coded
        config = self.get_config()
        self.assertEqual(config.server_config.data_locator__s3__region_name, "us-east-1")
        # incorrectly formatted
        datapath = "s3://shouldnt/work"
        file_name = self.custom_app_config(
            dataset_datapath=datapath, config_file_name=self.config_file_name, data_locater_region_name="true"
        )
        config = AppConfig()
        config.update_from_config_file(file_name)
        with self.assertRaises(ConfigurationError):
            config.server_config.handle_data_locator()

    def test_handle_app___can_use_envar_port(self):
        config = self.get_config(port=24)
        self.assertEqual(config.server_config.app__port, 24)

        # Note if the port is set in the config file it will NOT be overwritten by a different envvar
        os.environ["CXG_SERVER_PORT"] = "4008"
        self.config = AppConfig()
        self.config.update_server_config(app__flask_secret_key="secret")
        self.config.server_config.handle_app(self.context)
        self.assertEqual(self.config.server_config.app__port, 4008)
        del os.environ["CXG_SERVER_PORT"]

    def test_handle_app__can_get_secret_key_from_envvar_or_config_file_with_envvar_given_preference(self):
        config = self.get_config(flask_secret_key="KEY_FROM_FILE")
        self.assertEqual(config.server_config.app__flask_secret_key, "KEY_FROM_FILE")

        os.environ["CXG_SECRET_KEY"] = "KEY_FROM_ENV"
        config.external_config.handle_environment(self.context)
        self.assertEqual(config.server_config.app__flask_secret_key, "KEY_FROM_ENV")

    def test_config_for_single_dataset(self):
        file_name = self.custom_app_config(config_file_name="single_dataset.yml", dataset_datapath=f"{H5AD_FIXTURE}")
        config = AppConfig()
        config.update_from_config_file(file_name)
        config.server_config.handle_single_dataset(self.context)

        file_name = self.custom_app_config(
            config_file_name="single_dataset_with_about.yml",
            about="www.cziscience.com",
            dataset_datapath=f"{H5AD_FIXTURE}",
        )
        config = AppConfig()
        config.update_from_config_file(file_name)
        with self.assertRaises(ConfigurationError):
            config.server_config.handle_single_dataset(self.context)
