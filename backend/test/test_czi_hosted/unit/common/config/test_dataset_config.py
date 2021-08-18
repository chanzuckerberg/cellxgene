import json
import os
import tempfile

import unittest
from time import sleep
from unittest.mock import patch

from backend.czi_hosted.common.annotations.hosted_tiledb import AnnotationsHostedTileDB
from backend.czi_hosted.common.annotations.local_file_csv import AnnotationsLocalFile
from backend.czi_hosted.common.config.app_config import AppConfig
from backend.czi_hosted.common.config.base_config import BaseConfig
from backend.test import PROJECT_ROOT, FIXTURES_ROOT

from backend.common.errors import ConfigurationError
from backend.test.test_czi_hosted.unit.common.config import ConfigTests


class TestDatasetConfig(ConfigTests):
    def setUp(self):
        self.config_file_name = f"{unittest.TestCase.id(self).split('.')[-1]}.yml"
        self.config = AppConfig()
        self.config.update_server_config(app__flask_secret_key="secret")
        self.config.update_server_config(multi_dataset__dataroot=FIXTURES_ROOT)
        self.dataset_config = self.config.default_dataset_config
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

    def test_init_datatset_config_sets_vars_from_default_config(self):
        config = AppConfig()
        self.assertEqual(config.default_dataset_config.presentation__max_categories, 1000)
        self.assertEqual(config.default_dataset_config.user_annotations__type, "local_file_csv")
        self.assertEqual(config.default_dataset_config.diffexp__lfc_cutoff, 0.01)

    @patch("backend.czi_hosted.common.config.dataset_config.BaseConfig.validate_correct_type_of_configuration_attribute")
    def test_complete_config_checks_all_attr(self, mock_check_attrs):
        mock_check_attrs.side_effect = BaseConfig.validate_correct_type_of_configuration_attribute()
        self.dataset_config.complete_config(self.context)
        self.assertEqual(mock_check_attrs.call_count, 19)

    def test_app_sets_script_vars(self):
        config = self.get_config(scripts=["path/to/script"])
        config.default_dataset_config.handle_app()

        self.assertEqual(config.default_dataset_config.app__scripts, [{"src": "path/to/script"}])

        config = self.get_config(scripts=[{"src": "path/to/script", "more": "different/script/path"}])
        config.default_dataset_config.handle_app()
        self.assertEqual(
            config.default_dataset_config.app__scripts, [{"src": "path/to/script", "more": "different/script/path"}]
        )

        config = self.get_config(scripts=["path/to/script", "different/script/path"])
        config.default_dataset_config.handle_app()
        # TODO @madison -- is this the desired functionality?
        self.assertEqual(
            config.default_dataset_config.app__scripts, [{"src": "path/to/script"}, {"src": "different/script/path"}]
        )

        config = self.get_config(scripts=[{"more": "different/script/path"}])
        with self.assertRaises(ConfigurationError):
            config.default_dataset_config.handle_app()

    def test_handle_user_annotations_ensures_auth_is_enabled_with_valid_auth_type(self):
        config = self.get_config(enable_users_annotations="true", authentication_enable="false")
        config.server_config.complete_config(self.context)
        with self.assertRaises(ConfigurationError):
            config.default_dataset_config.handle_user_annotations(self.context)

        config = self.get_config(enable_users_annotations="true", authentication_enable="true", auth_type="pretend")
        with self.assertRaises(ConfigurationError):
            config.server_config.complete_config(self.context)

    def test_handle_user_annotations__adds_warning_message_if_annotation_vars_set_when_annotations_disabled(self):
        config = self.get_config(
            enable_users_annotations="false", authentication_enable="false", db_uri="shouldnt/be/set"
        )
        config.default_dataset_config.handle_user_annotations(self.context)

        self.assertEqual(self.context["messages"], ["Warning: db_uri ignored as annotations are disabled."])

    @patch("backend.czi_hosted.common.config.dataset_config.DbUtils")
    def test_handle_user_annotations__instantiates_user_annotations_class_correctly(self, mock_db_utils):
        mock_db_utils.return_value = "123"
        config = self.get_config(
            enable_users_annotations="true", authentication_enable="true", annotation_type="local_file_csv"
        )
        config.server_config.complete_config(self.context)
        config.default_dataset_config.handle_user_annotations(self.context)
        self.assertIsInstance(config.default_dataset_config.user_annotations, AnnotationsLocalFile)

        config = self.get_config(
            enable_users_annotations="true",
            authentication_enable="true",
            annotation_type="hosted_tiledb_array",
            db_uri="gotta/set/this",
            hosted_file_directory="and/this",
        )
        config.server_config.complete_config(self.context)
        config.default_dataset_config.handle_user_annotations(self.context)
        self.assertIsInstance(config.default_dataset_config.user_annotations, AnnotationsHostedTileDB)

        config = self.get_config(
            enable_users_annotations="true", authentication_enable="true", annotation_type="NOT_REAL"
        )
        config.server_config.complete_config(self.context)
        with self.assertRaises(ConfigurationError):
            config.default_dataset_config.handle_user_annotations(self.context)

    def test_handle_local_file_csv_annotations__sets_dir_if_not_passed_in(self):
        config = self.get_config(
            enable_users_annotations="true", authentication_enable="true", annotation_type="local_file_csv"
        )
        config.server_config.complete_config(self.context)
        config.default_dataset_config.handle_local_file_csv_annotations()
        self.assertIsInstance(config.default_dataset_config.user_annotations, AnnotationsLocalFile)
        cwd = os.getcwd()
        self.assertEqual(config.default_dataset_config.user_annotations._get_output_dir(), cwd)

    def test_handle_diffexp__raises_warning_for_large_datasets(self):
        config = self.get_config(lfc_cutoff=0.02, enable_difexp="true", top_n=15)
        config.server_config.complete_config(self.context)
        config.default_dataset_config.handle_diffexp(self.context)
        self.assertEqual(len(self.context["messages"]), 0)

    def test_multi_dataset(self):
        config = AppConfig()
        # test for illegal url_dataroots
        for illegal in ("../b", "!$*", "\\n", "", "(bad)"):
            config.update_server_config(
                app__flask_secret_key="secret",
                multi_dataset__dataroot={"tag": {"base_url": illegal, "dataroot": f"{PROJECT_ROOT}/example-dataset"}},
            )
            with self.assertRaises(ConfigurationError):
                config.complete_config()

        # test for legal url_dataroots
        for legal in ("d", "this.is-okay_", "a/b"):
            config.update_server_config(
                app__flask_secret_key="secret",
                multi_dataset__dataroot={"tag": {"base_url": legal, "dataroot": f"{PROJECT_ROOT}/example-dataset"}},
            )
            config.complete_config()

        # test that multi dataroots work end to end
        config.update_server_config(
            app__flask_secret_key="secret",
            multi_dataset__dataroot=dict(
                s1=dict(dataroot=f"{PROJECT_ROOT}/example-dataset", base_url="set1/1/2"),
                s2=dict(dataroot=f"{FIXTURES_ROOT}", base_url="set2"),
                s3=dict(dataroot=f"{FIXTURES_ROOT}", base_url="set3"),
            ),
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

        server = self.create_app(config)

        server.testing = True
        session = server.test_client()

        with self.subTest("Test config for dataroot /set1/1/2/ returns the s1 config"):
            response1 = session.get("/set1/1/2/pbmc3k.h5ad/api/v0.2/config")
            data_config_set_1 = json.loads(response1.data)

            self.assertEqual(data_config_set_1["config"]["displayNames"]["dataset"], "pbmc3k")
            self.assertFalse(data_config_set_1["config"]["parameters"]["annotations"])
            self.assertFalse(data_config_set_1["config"]["parameters"]["disable-diffexp"])
            self.assertEqual(data_config_set_1["config"]["parameters"]["about_legal_tos"], "tos_set1.html")

        with self.subTest("Test config for dataroot /set2 returns the s2 config"):
            response2 = session.get("/set2/pbmc3k.cxg/api/v0.2/config")
            data_config_set_2 = json.loads(response2.data)
            self.assertEqual(data_config_set_2["config"]["displayNames"]["dataset"], "pbmc3k")
            self.assertTrue(data_config_set_2["config"]["parameters"]["annotations"])
            self.assertEqual(data_config_set_2["config"]["parameters"]["about_legal_tos"], "tos_set2.html")

        with self.subTest("Test config for dataroot /set3/ returns the default dataset config"):
            response3 = session.get("/set3/pbmc3k.cxg/api/v0.2/config")
            data_config_set_3 = json.loads(response3.data)
            self.assertEqual(data_config_set_3["config"]["displayNames"]["dataset"], "pbmc3k")
            self.assertTrue(data_config_set_3["config"]["parameters"]["annotations"])
            self.assertFalse(data_config_set_3["config"]["parameters"]["disable-diffexp"])
            self.assertEqual(data_config_set_3["config"]["parameters"]["about_legal_tos"], "tos_default.html")


        response = session.get("/health")
        self.assertEqual(json.loads(response.data)["status"], "pass")

    def test_configfile_with_specialization(self):
        # test that per_dataset_config config load the default config, then the specialized config

        with tempfile.TemporaryDirectory() as tempdir:
            configfile = os.path.join(tempdir, "config.yaml")
            with open(configfile, "w") as fconfig:
                config = """
                server:
                    multi_dataset:
                        dataroot:
                            test:
                                base_url: test
                                dataroot: fake_dataroot

                dataset:
                    user_annotations:
                        enable: false
                        type: hosted_tiledb_array
                        hosted_tiledb_array:
                            db_uri: fake_db_uri
                            hosted_file_directory: fake_dir

                per_dataset_config:
                    test:
                        user_annotations:
                            enable: true
                """
                fconfig.write(config)

            app_config = AppConfig()
            app_config.update_from_config_file(configfile)

            test_config = app_config.dataroot_config["test"]

            # test config from default
            self.assertEqual(test_config.user_annotations__type, "hosted_tiledb_array")
            self.assertEqual(test_config.user_annotations__hosted_tiledb_array__db_uri, "fake_db_uri")

            # test config from specialization
            self.assertTrue(test_config.user_annotations__enable)
