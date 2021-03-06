import os
import tempfile
import unittest

import yaml

from local_server.default_config import default_config
from local_server.common.config.app_config import AppConfig
from local_server.test.unit.common.config import ConfigTests
from local_server.common.errors import ConfigurationError
from local_server.test import FIXTURES_ROOT, H5AD_FIXTURE


class AppConfigTest(ConfigTests):
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
            dataset_datapath=H5AD_FIXTURE, config_file_name=self.config_file_name, **kwargs
        )
        config = AppConfig()
        config.update_from_config_file(file_name)
        return config

    def test_get_default_config_correctly_reads_default_config_file(self):
        app_default_config = AppConfig().default_config

        expected_config = yaml.load(default_config, Loader=yaml.Loader)

        server_config = app_default_config["server"]
        dataset_config = app_default_config["dataset"]

        expected_server_config = expected_config["server"]
        expected_dataset_config = expected_config["dataset"]

        self.assertDictEqual(app_default_config, expected_config)
        self.assertDictEqual(server_config, expected_server_config)
        self.assertDictEqual(dataset_config, expected_dataset_config)

    def test_get_dataset_config_returns_dataset_config_for_single_datasets(self):
        datapath = f"{FIXTURES_ROOT}/1e4dfec4-c0b2-46ad-a04e-ff3ffb3c0a8f.h5ad"
        file_name = self.custom_app_config(dataset_datapath=datapath, config_file_name=self.config_file_name)
        config = AppConfig()
        config.update_from_config_file(file_name)

        self.assertEqual(config.get_dataset_config(), config.dataset_config)

    def test_update_server_config_updates_server_config_and_config_status(self):
        config = self.get_config()
        config.complete_config()
        config.check_config()
        config.update_server_config(single_dataset__datapath=H5AD_FIXTURE)
        with self.assertRaises(ConfigurationError):
            config.server_config.check_config()

    def test_write_config_outputs_yaml_with_all_config_vars(self):
        config = self.get_config()
        config.write_config(f"{FIXTURES_ROOT}/tmp_dir/write_config.yml")
        with open(f"{FIXTURES_ROOT}/tmp_dir/{self.config_file_name}", "r") as default_config:
            default_config_yml = yaml.safe_load(default_config)

        with open(f"{FIXTURES_ROOT}/tmp_dir/write_config.yml", "r") as output_config:
            output_config_yml = yaml.safe_load(output_config)
        self.maxDiff = None
        self.assertEqual(default_config_yml, output_config_yml)

    def test_update_app_config(self):
        config = AppConfig()
        config.update_server_config(app__verbose=True, single_dataset__datapath="datapath")
        vars = config.server_config.changes_from_default()
        self.assertCountEqual(vars, [("app__verbose", True, False), ("single_dataset__datapath", "datapath", None)])

        config = AppConfig()
        config.update_dataset_config(app__scripts=(), app__inline_scripts=())
        vars = config.server_config.changes_from_default()
        self.assertCountEqual(vars, [])

        config = AppConfig()
        config.update_dataset_config(app__scripts=[], app__inline_scripts=[])
        vars = config.dataset_config.changes_from_default()
        self.assertCountEqual(vars, [])

        config = AppConfig()
        config.update_dataset_config(app__scripts=("a", "b"), app__inline_scripts=["c", "d"])
        vars = config.dataset_config.changes_from_default()
        self.assertCountEqual(vars, [("app__scripts", ["a", "b"], []), ("app__inline_scripts", ["c", "d"], [])])

    def test_configfile_no_server_section(self):
        # test a config file without a dataset section

        with tempfile.TemporaryDirectory() as tempdir:
            configfile = os.path.join(tempdir, "config.yaml")
            with open(configfile, "w") as fconfig:
                config = """
                dataset:
                    user_annotations:
                        enable: false
                """
                fconfig.write(config)

            app_config = AppConfig()
            app_config.update_from_config_file(configfile)
            server_changes = app_config.server_config.changes_from_default()
            dataset_changes = app_config.dataset_config.changes_from_default()
            self.assertEqual(server_changes, [])
            self.assertEqual(dataset_changes, [("user_annotations__enable", False, True)])

    def test_simple_update_single_config_from_path_and_value(self):
        """Update a simple config parameter"""

        config = AppConfig()
        config.server_config.single_dataset__datapath = "my/data/path"

        # test simple value in server
        config.update_single_config_from_path_and_value(["server", "app", "flask_secret_key"], "mysecret")
        self.assertEqual(config.server_config.app__flask_secret_key, "mysecret")

        # test simple value in default dataset
        config.update_single_config_from_path_and_value(
            ["dataset", "user_annotations", "ontology", "obo_location"], "dummy_location",
        )
        self.assertEqual(config.dataset_config.user_annotations__ontology__obo_location, "dummy_location")

        # error checking
        bad_paths = [
            (
                ["dataset", "does", "not", "exist"],
                "unknown config parameter at path: '['dataset', 'does', 'not', 'exist']'",
            ),
            (["does", "not", "exist"], "path must start with 'server', or 'dataset'"),
            ([], "path must start with 'server', or 'dataset'"),
            ([1, 2, 3], "path must be a list of strings, got '[1, 2, 3]'"),
            ("string", "path must be a list of strings, got 'string'"),
        ]
        for bad_path, error_message in bad_paths:
            with self.assertRaises(ConfigurationError) as config_error:
                config.update_single_config_from_path_and_value(bad_path, "value")

            self.assertEqual(config_error.exception.message, error_message)
