import os
import unittest
from unittest import mock
from unittest.mock import patch

from local_server.common.config.base_config import BaseConfig
from local_server.common.utils.utils import find_available_port
from local_server.test import PROJECT_ROOT, FIXTURES_ROOT

import requests

from local_server.common.config.app_config import AppConfig
from local_server.common.errors import ConfigurationError
from local_server.test import test_server
from local_server.test.unit.common.config import ConfigTests


def mockenv(**envvars):
    return mock.patch.dict(os.environ, envvars)


class TestServerConfig(ConfigTests):
    def setUp(self):
        self.config_file_name = f"{unittest.TestCase.id(self).split('.')[-1]}.yml"
        self.config = AppConfig()
        self.config.update_server_config(app__flask_secret_key="secret")
        self.config.update_server_config(multi_dataset__dataroot=FIXTURES_ROOT)
        self.server_config = self.config.server_config
        self.config.complete_config()

        message_list = []

        def noop(message):
            message_list.append(message)

        messagefn = noop
        self.context = dict(messagefn=messagefn, messages=message_list)

    def get_config(self, **kwargs):
        file_name = self.custom_app_config(
            dataroot=f"{FIXTURES_ROOT}", config_file_name=self.config_file_name, **kwargs
        )
        config = AppConfig()
        config.update_from_config_file(file_name)
        return config

    def test_init_raises_error_if_default_config_is_invalid(self):
        invalid_config = self.get_config(port="not_valid")
        with self.assertRaises(ConfigurationError):
            invalid_config.complete_config()

    @patch("local_server.common.config.server_config.BaseConfig.validate_correct_type_of_configuration_attribute")
    def test_complete_config_checks_all_attr(self, mock_check_attrs):
        mock_check_attrs.side_effect = BaseConfig.validate_correct_type_of_configuration_attribute()
        self.server_config.complete_config(self.context)
        self.assertEqual(mock_check_attrs.call_count, 30)

    def test_handle_app__throws_error_if_port_doesnt_exist(self):
        config = self.get_config(port=99999999)
        with self.assertRaises(ConfigurationError):
            config.server_config.handle_app(self.context)

    @patch("local_server.common.config.server_config.discover_s3_region_name")
    def test_handle_data_locator_works_for_default_types(self, mock_discover_region_name):
        mock_discover_region_name.return_value = None
        # Default config
        self.assertEqual(self.config.server_config.data_locator__s3__region_name, None)
        # hard coded
        config = self.get_config()
        self.assertEqual(config.server_config.data_locator__s3__region_name, "us-east-1")
        # incorrectly formatted
        dataroot = {
            "d1": {"base_url": "set1", "dataroot": "/path/to/set1_datasets/"},
            "d2": {"base_url": "set2/subdir", "dataroot": "s3://shouldnt/work"},
        }
        file_name = self.custom_app_config(
            dataroot=dataroot, config_file_name=self.config_file_name, data_locater_region_name="true"
        )
        config = AppConfig()
        config.update_from_config_file(file_name)
        with self.assertRaises(ConfigurationError):
            config.server_config.handle_data_locator()

    @patch("local_server.common.config.server_config.discover_s3_region_name")
    def test_handle_data_locator_can_read_from_dataroot(self, mock_discover_region_name):
        mock_discover_region_name.return_value = "us-west-2"
        dataroot = {
            "d1": {"base_url": "set1", "dataroot": "/path/to/set1_datasets/"},
            "d2": {"base_url": "set2/subdir", "dataroot": "s3://hosted-cellxgene-dev"},
        }
        file_name = self.custom_app_config(
            dataroot=dataroot, config_file_name=self.config_file_name, data_locater_region_name="true"
        )
        config = AppConfig()
        config.update_from_config_file(file_name)
        config.server_config.handle_data_locator()
        self.assertEqual(config.server_config.data_locator__s3__region_name, "us-west-2")
        mock_discover_region_name.assert_called_once_with("s3://hosted-cellxgene-dev")

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

    def test_handle_app__sets_web_base_url(self):
        config = self.get_config(web_base_url="anything.com")
        self.assertEqual(config.server_config.app__web_base_url, "anything.com")

    def test_handle_data_source__errors_when_passed_zero_or_two_dataroots(self):
        file_name = self.custom_app_config(
            dataroot=f"{FIXTURES_ROOT}",
            config_file_name="two_data_roots.yml",
            dataset_datapath=f"{FIXTURES_ROOT}/pbmc3k-CSC-gz.h5ad",
        )
        config = AppConfig()
        config.update_from_config_file(file_name)
        with self.assertRaises(ConfigurationError):
            config.server_config.handle_data_source()

        file_name = self.custom_app_config(config_file_name="zero_roots.yml")
        config = AppConfig()
        config.update_from_config_file(file_name)
        with self.assertRaises(ConfigurationError):
            config.server_config.handle_data_source()

    def test_get_api_base_url_works(self):

        # test the api_base_url feature, and that it can contain a path
        config = AppConfig()
        backend_port = find_available_port("localhost", 10000)
        config.update_server_config(
            app__flask_secret_key="secret",
            app__api_base_url=f"http://localhost:{backend_port}/additional/path",
            multi_dataset__dataroot=f"{PROJECT_ROOT}/example-dataset",
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

    def test_get_web_base_url_works(self):
        config = self.get_config(web_base_url="www.thisisawebsite.com")
        web_base_url = config.server_config.get_web_base_url()
        self.assertEqual(web_base_url, "www.thisisawebsite.com")

        config = self.get_config(web_base_url="local", port=12)
        web_base_url = config.server_config.get_web_base_url()
        self.assertEqual(web_base_url, "http://localhost:12")

        config = self.get_config(web_base_url="www.thisisawebsite.com/")
        web_base_url = config.server_config.get_web_base_url()
        self.assertEqual(web_base_url, "www.thisisawebsite.com")

        config = self.get_config(api_base_url="www.api_base.com/")
        web_base_url = config.server_config.get_web_base_url()
        self.assertEqual(web_base_url, "www.api_base.com")

    def test_config_for_single_dataset(self):
        file_name = self.custom_app_config(
            config_file_name="single_dataset.yml", dataset_datapath=f"{FIXTURES_ROOT}/pbmc3k.cxg"
        )
        config = AppConfig()
        config.update_from_config_file(file_name)
        config.server_config.handle_single_dataset(self.context)
        self.assertIsNotNone(config.server_config.matrix_data_cache_manager)

        file_name = self.custom_app_config(
            config_file_name="single_dataset_with_about.yml",
            about="www.cziscience.com",
            dataset_datapath=f"{FIXTURES_ROOT}/pbmc3k.cxg",
        )
        config = AppConfig()
        config.update_from_config_file(file_name)
        with self.assertRaises(ConfigurationError):
            config.server_config.handle_single_dataset(self.context)

    def test_multi_dataset_raises_error_for_illegal_routes(self):
        # test for illegal url_dataroots
        for illegal in ("../b", "!$*", "\\n", "", "(bad)"):
            self.config.update_server_config(
                multi_dataset__dataroot={"tag": {"base_url": illegal, "dataroot": f"{PROJECT_ROOT}/example-dataset"}}
            )
            with self.assertRaises(ConfigurationError):
                self.config.complete_config()

    def test_multidataset_works_for_legal_routes(self):
        # test for legal url_dataroots
        for legal in ("d", "this.is-okay_", "a/b"):
            self.config.update_server_config(
                multi_dataset__dataroot={"tag": {"base_url": legal, "dataroot": f"{PROJECT_ROOT}/example-dataset"}}
            )
            self.config.complete_config()

    def test_mulitdatasets_work_e2e(self):
        # test that multi dataroots work end to end
        self.config.update_server_config(
            multi_dataset__dataroot=dict(
                s1=dict(dataroot=f"{PROJECT_ROOT}/example-dataset", base_url="set1/1/2"),
                s2=dict(dataroot=f"{FIXTURES_ROOT}", base_url="set2"),
                s3=dict(dataroot=f"{FIXTURES_ROOT}", base_url="set3"),
            )
        )

        # Change this default to test if the dataroot overrides below work.
        self.config.update_default_dataset_config(app__about_legal_tos="tos_default.html")

        # specialize the configs for set1
        self.config.add_dataroot_config(
            "s1", user_annotations__enable=False, diffexp__enable=True, app__about_legal_tos="tos_set1.html"
        )

        # specialize the configs for set2
        self.config.add_dataroot_config(
            "s2", user_annotations__enable=True, diffexp__enable=False, app__about_legal_tos="tos_set2.html"
        )

        # no specializations for set3 (they get the default dataset config)
        self.config.complete_config()

        with test_server(app_config=self.config) as server:
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

            # access a dataset (no slash)
            response = session.get(f"{server}/set2/pbmc3k.cxg")
            self.assertEqual(response.status_code, 200)

            # access a dataset (with slash)
            response = session.get(f"{server}/set2/pbmc3k.cxg/")
            self.assertEqual(response.status_code, 200)

    @patch("local_server.data_cxg.cxg_adaptor.CxgAdaptor.set_tiledb_context")
    def test_handle_adaptor(self, mock_tiledb_context):
        custom_config = self.custom_app_config(
            dataroot=f"{FIXTURES_ROOT}", cxg_tile_cache_size=10, cxg_num_reader_threads=2
        )
        config = AppConfig()
        config.update_from_config_file(custom_config)
        config.server_config.handle_adaptor()
        mock_tiledb_context.assert_called_once_with(
            {"sm.tile_cache_size": 10, "sm.num_reader_threads": 2, "vfs.s3.region": "us-east-1"}
        )
