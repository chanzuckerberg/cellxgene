import unittest
import json

from server.common.utils.data_locator import DataLocator
from server.data_soma.soma_adaptor import SomaAdaptor
from server.common.config.app_config import AppConfig
from test import PROJECT_ROOT


class SomaDataLoadAdaptorTest(unittest.TestCase):
    """
    Test file loading, including deferred loading/update.
    """

    def setUp(self):
        self.data_file = DataLocator(f"{PROJECT_ROOT}/example-dataset/tiledb-data/pbmc3k_processed")
        config = AppConfig()
        config.update_server_config(single_dataset__datapath=self.data_file.path)
        config.update_server_config(app__flask_secret_key="secret")
        config.complete_config()
        self.data = SomaAdaptor(self.data_file, config)

    def test_delayed_load_data(self):
        self.data._create_schema()
        self.assertEqual(self.data.cell_count, 2638)
        self.assertEqual(self.data.gene_count, 1838)
        epsilon = 0.000_005
        self.assertTrue(self.data.data.X.data.df().to_numpy()[0, 0] - -0.171_469_51 < epsilon)

    def test_diffexp_topN(self):
        f1 = {"filter": {"obs": {"index": [[0, 500]]}}}
        f2 = {"filter": {"obs": {"index": [[500, 1000]]}}}
        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"]))
        self.assertEqual(len(result["positive"]), 10)
        self.assertEqual(len(result["negative"]), 10)
        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"], 20))
        self.assertEqual(len(result["positive"]), 20)
        self.assertEqual(len(result["negative"]), 20)


class SomaDataLocatorAdaptorTest(unittest.TestCase):
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
        """run these each time we load the data"""
        self.assertIsNotNone(data)
        self.assertEqual(data.cell_count, 2638)
        self.assertEqual(data.gene_count, 1838)

    def test_posix_file(self):
        locator = DataLocator("example-dataset/tiledb-data/pbmc3k_processed")
        self.assertEqual(SomaAdaptor.file_size(locator), 480)
        config = self.get_basic_config()
        config.update_server_config(single_dataset__datapath=locator.path)
        config.complete_config()
        data = SomaAdaptor(locator, config)
        self.stdAsserts(data)

    # TODO: using SOMA with a TileDB URI is not tested yet
    # def test_url_https(self):
    #     url = "https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad"
    #     locator = DataLocator(url)
    #     config = self.get_basic_config()
    #     data = AnndataAdaptor(locator, config)
    #     self.stdAsserts(data)

    # def test_url_http(self):
    #     url = "http://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad"
    #     locator = DataLocator(url)
    #     config = self.get_basic_config()
    #     data = AnndataAdaptor(locator, config)
    #     self.stdAsserts(data)
