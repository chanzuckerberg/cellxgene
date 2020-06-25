import unittest
from server.common.app_config import AppConfig
from server.common.errors import ConfigurationError
from server.test import PROJECT_ROOT, test_server
import requests

# NOTE, there are more tests that should be written for AppConfig.
# this is just a start.


class AppConfigTest(unittest.TestCase):
    def test_update(self):
        c = AppConfig()
        c.update_server_config(app__verbose=True, multi_dataset__dataroot="datadir")
        v = c.server_config.changes_from_default()
        self.assertCountEqual(v, [("app__verbose", True, False), ("multi_dataset__dataroot", "datadir", None)])

        c = AppConfig()
        c.update_default_dataset_config(app__scripts=(), app__inline_scripts=())
        v = c.server_config.changes_from_default()
        self.assertCountEqual(v, [])

        c = AppConfig()
        c.update_default_dataset_config(app__scripts=[], app__inline_scripts=[])
        v = c.default_dataset_config.changes_from_default()
        self.assertCountEqual(v, [])

        c = AppConfig()
        c.update_default_dataset_config(app__scripts=("a", "b"), app__inline_scripts=["c", "d"])
        v = c.default_dataset_config.changes_from_default()
        self.assertCountEqual(v, [("app__scripts", ["a", "b"], []), ("app__inline_scripts", ["c", "d"], [])])

    def test_multi_dataset(self):

        c = AppConfig()
        # test for illegal url_dataroots
        for illegal in ("a/b", "../b", "!$*", "\\n", "", "(bad)"):
            c.update_server_config(multi_dataset__dataroot={illegal: f"{PROJECT_ROOT}/example-dataset"})
            with self.assertRaises(ConfigurationError):
                c.complete_config()

        # test for legal url_dataroots
        for legal in (
            "d",
            "this.is-okay_",
        ):
            c.update_server_config(multi_dataset__dataroot={legal: f"{PROJECT_ROOT}/example-dataset"})
            c.complete_config()

        # test that multi dataroots work end to end
        c.update_server_config(
            multi_dataset__dataroot=dict(
                set1=f"{PROJECT_ROOT}/example-dataset", set2=f"{PROJECT_ROOT}/server/test/test_datasets"
            )
        )
        c.complete_config()

        with test_server(app_config=c) as server:
            session = requests.Session()

            r = session.get(f"{server}/set1/pbmc3k.h5ad/api/v0.2/config")
            data_config = r.json()
            assert data_config["config"]["displayNames"]["dataset"] == "pbmc3k"

            r = session.get(f"{server}/set2/pbmc3k.cxg/api/v0.2/config")
            data_config = r.json()
            assert data_config["config"]["displayNames"]["dataset"] == "pbmc3k"

            r = session.get(f"{server}/health")
            assert r.json()["status"] == "pass"
