import os
import tempfile

import unittest
from unittest.mock import patch

from backend.server.common.annotations.local_file_csv import AnnotationsLocalFile
from backend.server.common.config.app_config import AppConfig
from backend.server.common.config.base_config import BaseConfig
from backend.test import FIXTURES_ROOT, H5AD_FIXTURE

from backend.common.errors import ConfigurationError
from backend.test.test_server.unit.common.config import ConfigTests


class TestDatasetConfig(ConfigTests):
    def setUp(self):
        self.config_file_name = f"{unittest.TestCase.id(self).split('.')[-1]}.yml"
        self.config = AppConfig()
        self.config.update_server_config(app__flask_secret_key="secret")
        self.config.update_server_config(single_dataset__datapath=H5AD_FIXTURE)
        self.dataset_config = self.config.dataset_config
        self.config.complete_config()
        message_list = []

        def noop(message):
            message_list.append(message)

        messagefn = noop
        self.context = dict(messagefn=messagefn, messages=message_list)

    def get_config(self, **kwargs):
        file_name = self.custom_app_config(dataset_datapath=H5AD_FIXTURE, **kwargs)
        config = AppConfig()
        config.update_from_config_file(file_name)
        return config

    def test_init_datatset_config_sets_vars_from_config(self):
        config = AppConfig()
        self.assertEqual(config.dataset_config.presentation__max_categories, 1000)
        self.assertEqual(config.dataset_config.user_annotations__type, "local_file_csv")
        self.assertEqual(config.dataset_config.diffexp__lfc_cutoff, 0.01)

    @patch("backend.server.common.config.dataset_config.BaseConfig.validate_correct_type_of_configuration_attribute")
    def test_complete_config_checks_all_attr(self, mock_check_attrs):
        mock_check_attrs.side_effect = BaseConfig.validate_correct_type_of_configuration_attribute()
        self.dataset_config.complete_config(self.context)
        self.assertIsNotNone(self.config.server_config.data_adaptor)
        self.assertEqual(mock_check_attrs.call_count, 16)

    def test_app_sets_script_vars(self):
        config = self.get_config(scripts=["path/to/script"])
        config.dataset_config.handle_app()

        self.assertEqual(config.dataset_config.app__scripts, [{"src": "path/to/script"}])

        config = self.get_config(scripts=[{"src": "path/to/script", "more": "different/script/path"}])
        config.dataset_config.handle_app()
        self.assertEqual(
            config.dataset_config.app__scripts, [{"src": "path/to/script", "more": "different/script/path"}]
        )

        config = self.get_config(scripts=["path/to/script", "different/script/path"])
        config.dataset_config.handle_app()
        # TODO @madison -- is this the desired functionality?
        self.assertEqual(
            config.dataset_config.app__scripts, [{"src": "path/to/script"}, {"src": "different/script/path"}]
        )

        config = self.get_config(scripts=[{"more": "different/script/path"}])
        with self.assertRaises(ConfigurationError):
            config.dataset_config.handle_app()

    def test_handle_user_annotations_ensures_auth_is_enabled_with_valid_auth_type(self):
        config = self.get_config(enable_users_annotations="true", authentication_enable="false")
        config.server_config.complete_config(self.context)
        with self.assertRaises(ConfigurationError):
            config.dataset_config.handle_user_annotations(self.context)

        config = self.get_config(enable_users_annotations="true", authentication_enable="true", auth_type="pretend")
        with self.assertRaises(ConfigurationError):
            config.server_config.complete_config(self.context)

    def test_handle_user_annotations__instantiates_user_annotations_class_correctly(self):
        config = self.get_config(
            enable_users_annotations="true", authentication_enable="true", annotation_type="local_file_csv"
        )
        config.server_config.complete_config(self.context)
        config.dataset_config.handle_user_annotations(self.context)
        self.assertIsInstance(config.dataset_config.user_annotations, AnnotationsLocalFile)

        config = self.get_config(
            enable_users_annotations="true", authentication_enable="true", annotation_type="NOT_REAL"
        )
        config.server_config.complete_config(self.context)
        with self.assertRaises(ConfigurationError):
            config.dataset_config.handle_user_annotations(self.context)

    def test_handle_local_file_csv_annotations__sets_dir_if_not_passed_in(self):
        config = self.get_config(
            enable_users_annotations="true", authentication_enable="true", annotation_type="local_file_csv"
        )
        config.server_config.complete_config(self.context)
        config.dataset_config.handle_local_file_csv_annotations(self.context)
        self.assertIsInstance(config.dataset_config.user_annotations, AnnotationsLocalFile)
        cwd = os.getcwd()
        self.assertEqual(config.dataset_config.user_annotations._get_output_dir(), cwd)

    def test_handle_diffexp__raises_warning_for_large_datasets(self):
        config = self.get_config(lfc_cutoff=0.02, enable_difexp="true", top_n=15)
        config.server_config.complete_config(self.context)
        config.dataset_config.handle_diffexp(self.context)
        self.assertEqual(len(self.context["messages"]), 1)

    def test_configfile_with_specialization(self):
        # test that per_dataset_config config load the default config, then the specialized config

        with tempfile.TemporaryDirectory() as tempdir:
            configfile = os.path.join(tempdir, "config.yaml")
            with open(configfile, "w") as fconfig:
                config = """
                server:
                    single_dataset:
                        datapath: fake_datapath
                dataset:
                    user_annotations:
                        enable: false
                        type: local_file_csv
                        local_file_csv:
                            file: fake_file
                            directory: fake_dir
                """
                fconfig.write(config)

            app_config = AppConfig()
            app_config.update_from_config_file(configfile)

            test_config = app_config.dataset_config

            # test config from default
            self.assertEqual(test_config.user_annotations__type, "local_file_csv")
            self.assertEqual(test_config.user_annotations__local_file_csv__file, "fake_file")
