import os
import tempfile
import unittest

from server.common.config.app_config import AppConfig


class AppConfigTest(unittest.TestCase):

    def test_get_default_config_successfully_reads_yaml_file(self):
        pass

    def test_get_default_config_raises_an_error_for_invalid_config_files(self):
        pass

    def test_init_app_config_pulls_in_default_server_and_dataset_configurations(self):
        pass

    def test_get_dataset_config_handles_single_datasets(self):
        pass

    def test_get_dataset_config__returns_default_if_not_set(self):
        pass

    def test_get_dataset_config__returns_passed_in_dataset_config(self):
        pass

    def test_check_config_verifies_attributes_in_all_configs_are_checked(self):
        pass

    def test_update_server_config_updates_server_config_and_config_status(self):
        pass

    # Todo -- is this the desired/expected functionality?
    def test_update_default_datasets_updates_for_all_dataroots(self):
        pass

    def test_update_from_config_file_correctly_updates_server_and_datasets_configs(self):
        pass

    def test_write_config_outputs_yaml_with_all_config_vars(self):
        pass

    # Todo -- should this return non default dataset configurations for additional dataroots?
    def test_changes_from_default_returns_attributes_from_server_and_dataset_that_are_non_default(self):
        # Note this wont check additional datasets bc there are no default values for those
        pass

    def test_add_dataroot_config_creates_config_based_on_dataset_default_and_passed_in_params(self):
        pass

    def test_complete_config_checks_server_dataset_and_all_dataroot_configs(self):
        pass

    def test_update_app_config(self):
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

    def test_configfile_no_dataset_section(self):
        # test a config file without a dataset section

        with tempfile.TemporaryDirectory() as tempdir:
            configfile = os.path.join(tempdir, "config.yaml")
            with open(configfile, "w") as fconfig:
                config = """
                server:
                    multi_dataset:
                        dataroot: test_dataroot

                """
                fconfig.write(config)

            app_config = AppConfig()
            app_config.update_from_config_file(configfile)
            server_changes = app_config.server_config.changes_from_default()
            dataset_changes = app_config.default_dataset_config.changes_from_default()
            self.assertEqual(server_changes, [("multi_dataset__dataroot", "test_dataroot", None)])
            self.assertEqual(dataset_changes, [])

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
            dataset_changes = app_config.default_dataset_config.changes_from_default()
            self.assertEqual(server_changes, [])
            self.assertEqual(dataset_changes, [("user_annotations__enable", False, True)])

