import unittest
import json

from backend.common.utils.data_locator import DataLocator
from backend.server.data_anndata.anndata_adaptor import AnndataAdaptor
from backend.server.common.config.app_config import AppConfig
from backend.test import PROJECT_ROOT


class DataLoadAdaptorTest(unittest.TestCase):
    """
    Test file loading, including deferred loading/update.
    """

    def setUp(self):
        self.data_file = DataLocator(f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad")
        config = AppConfig()
        config.update_server_config(single_dataset__datapath=self.data_file.path)
        config.update_server_config(app__flask_secret_key="secret")
        config.complete_config()
        self.data = AnndataAdaptor(self.data_file, config)

    def test_delayed_load_data(self):
        self.data._create_schema()
        self.assertEqual(self.data.cell_count, 2638)
        self.assertEqual(self.data.gene_count, 1838)
        epsilon = 0.000_005
        self.assertTrue(self.data.data.X[0, 0] - -0.171_469_51 < epsilon)

    def test_diffexp_topN(self):
        f1 = {"filter": {"obs": {"index": [[0, 500]]}}}
        f2 = {"filter": {"obs": {"index": [[500, 1000]]}}}
        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"]))
        self.assertEqual(len(result["positive"]), 10)
        self.assertEqual(len(result["negative"]), 10)
        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"], 20))
        self.assertEqual(len(result["positive"]), 20)
        self.assertEqual(len(result["negative"]), 20)


class DataLocatorAdaptorTest(unittest.TestCase):
    """
    Test various types of data locators we expect to consume
    """

    def get_basic_config(self):
        config = AppConfig()
        config.update_server_config(
            single_dataset__obs_names=None,
            single_dataset__var_names=None,
        )
        config.update_server_config(app__flask_secret_key="secret")
        config.update_dataset_config(
            embeddings__names=["umap"],
            presentation__max_categories=100,
            diffexp__lfc_cutoff=0.01,
        )
        return config

    def stdAsserts(self, data):
        """ run these each time we load the data """
        self.assertIsNotNone(data)
        self.assertEqual(data.cell_count, 2638)
        self.assertEqual(data.gene_count, 1838)

    def test_posix_file(self):
        locator = DataLocator("../../example-dataset/pbmc3k.h5ad")
        config = self.get_basic_config()
        config.update_server_config(single_dataset__datapath=locator.path)
        config.complete_config()
        data = AnndataAdaptor(locator, config)
        self.stdAsserts(data)

    def test_url_https(self):
        url = "https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad"
        locator = DataLocator(url)
        config = self.get_basic_config()
        data = AnndataAdaptor(locator, config)
        self.stdAsserts(data)

    def test_url_http(self):
        url = "http://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad"
        locator = DataLocator(url)
        config = self.get_basic_config()
        data = AnndataAdaptor(locator, config)
        self.stdAsserts(data)
