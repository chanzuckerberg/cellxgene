import unittest

from local_server.common.config.app_config import AppConfig
from local_server.test import H5AD_FIXTURE
from local_server.test.unit.common.config import ConfigTests
from local_server.common.errors import ConfigurationError


class BaseConfigTest(ConfigTests):
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

    def test_mapping_creation_returns_map_of_server_and_dataset_config(self):
        config = AppConfig()
        mapping = config.dataset_config.create_mapping(config.default_config)
        self.assertIsNotNone(mapping["server__app__verbose"])
        self.assertIsNotNone(mapping["dataset__presentation__max_categories"])
        self.assertIsNotNone(mapping["dataset__user_annotations__ontology__obo_location"])

    def test_changes_from_default_returns_list_of_nondefault_config_values(self):
        config = self.get_config(verbose="true", lfc_cutoff=0.05)
        server_changes = config.server_config.changes_from_default()
        dataset_changes = config.dataset_config.changes_from_default()

        self.assertEqual(
            server_changes,
            [
                ("app__verbose", True, False),
                ("app__flask_secret_key", "secret", None),
                ("single_dataset__datapath", H5AD_FIXTURE, None),
                ('data_locator__s3__region_name', 'us-east-1', True)
            ],
        )
        self.assertEqual(dataset_changes, [("diffexp__lfc_cutoff", 0.05, 0.01)])

    def test_check_config_throws_error_if_attr_has_not_been_checked(self):
        config = self.get_config(verbose="true")
        config.complete_config()
        config.check_config()
        config.update_server_config(app__verbose=False)
        with self.assertRaises(ConfigurationError):
            config.check_config()
