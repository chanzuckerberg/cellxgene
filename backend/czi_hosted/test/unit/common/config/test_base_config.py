import unittest

from backend.czi_hosted.common.config.app_config import AppConfig
from backend.czi_hosted.test import FIXTURES_ROOT
from backend.common_utils.errors import ConfigurationError
from backend.czi_hosted.test.unit.common.config import ConfigTests


class BaseConfigTest(ConfigTests):
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

    def test_mapping_creation_returns_map_of_server_and_dataset_config(self):
        config = AppConfig()
        mapping = config.default_dataset_config.create_mapping(config.default_config)
        self.assertIsNotNone(mapping["server__app__verbose"])
        self.assertIsNotNone(mapping["dataset__presentation__max_categories"])
        self.assertIsNotNone(mapping["dataset__user_annotations__ontology__obo_location"])
        self.assertIsNotNone(mapping["server__multi_dataset__allowed_matrix_types"])

    def test_changes_from_default_returns_list_of_nondefault_config_values(self):
        config = self.get_config(verbose="true", lfc_cutoff=0.05)
        server_changes = config.server_config.changes_from_default()
        dataset_changes = config.default_dataset_config.changes_from_default()

        self.assertEqual(
            server_changes,
            [
                ("app__verbose", True, False),
                ("app__flask_secret_key", "secret", None),
                ("multi_dataset__dataroot", FIXTURES_ROOT, None),
                ("multi_dataset__matrix_cache__timelimit_s", 5, 30),
                ("data_locator__s3__region_name", "us-east-1", True),
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
