import os
import unittest
from unittest import mock
from unittest.mock import patch

import yaml

from server.common.config.base_config import BaseConfig
from server.common.utils.utils import find_available_port
from server.test import PROJECT_ROOT, test_server, FIXTURES_ROOT

import requests

from server.common.config.app_config import AppConfig
from server.common.errors import ConfigurationError
from server.test import test_server


def mockenv(**envvars):
    return mock.patch.dict(os.environ, envvars)


class TestServerConfig(unittest.TestCase):
    def setUp(self):
        self.config = AppConfig()
        # cls.default_config = get_default_config_for_testing()

        self.server_config = self.config.server_config
        self.config.update_server_config(
            multi_dataset__dataroot=dict(
                s2=dict(dataroot=f"{FIXTURES_ROOT}", base_url="set2"),
                s3=dict(dataroot=f"{FIXTURES_ROOT}", base_url="set3"),
            ),
            authentication__type="test",
        )

        def noop():
            pass
        self.config.complete_config()

        messagefn = noop
        self.context = dict(messagefn=messagefn)

    def test_init_raises_error_if_default_config_is_invalid(self):
        invalid_config = f"{FIXTURES_ROOT}/invalid_config.yml"
        with self.assertRaises(ConfigurationError):
            self.config.update_from_config_file(invalid_config)


    @patch('server.common.config.server_config.BaseConfig.check_attr')
    def test_complete_config_checks_all_attr(self, mock_check_attrs):
        mock_check_attrs.side_effect = BaseConfig.check_attr()
        self.server_config.complete_config(self.context)
        self.assertEqual(mock_check_attrs.call_count, 40)

    def test_handle_app__throws_error_if_port_in_use_or_doesnt_exist(self):
        pass

    def test_handle_app___can_use_envar_port(self):
        pass

    def test_handle_app__can_get_secret_key_from_envvar_or_config_file_with_envvar_given_preference(self):
        pass

    # todo ask seve about best way to test
    def test_csp_directives(self):
        pass

    def test_handle_app__sets_web_base_url(self):
        pass

    def test_handle_auth__gets_client_secret_from_envvars_or_config_with_envvars_given_preference(self):
        pass

    def test_exceeds_limit_returns_correctly(self):
        pass

    def test_handle_auth__correctly_creates_factories(self):
        pass

    def test_handle_data_locator__correctly_finds_s3_region(self):
        pass

    def test_handle_data_source__errors_when_passed_zero_or_two_dataroots(self):
        pass

    def test_get_api_base_url_works(self):

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

    def test_get_web_base_url_works(self):
        pass

    def test_config_for_single_dataset(self):
        pass

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

    def test_handle_diffexp(self):
        pass

    def test_handle_adaptor(self):
        pass
