import unittest
import json

from server.data_anndata.anndata_adaptor import AnndataAdaptor
from server.common.data_locator import DataLocator
from server.common.app_config import AppConfig
from server.test import PROJECT_ROOT


class DataLoadAdaptorTest(unittest.TestCase):
    """
    Test file loading, including deferred loading/update.
    """

    def setUp(self):
        self.data_file = DataLocator(f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad")
        config = AppConfig()
        config.update(single_dataset__datapath=self.data_file.path)
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
        self.assertEqual(len(result), 10)
        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"], 20))
        self.assertEqual(len(result), 20)


class DataLocatorAdaptorTest(unittest.TestCase):
    """
    Test various types of data locators we expect to consume
    """

    def setUp(self):
        self.args = {
            "embeddings__names": ["umap"],
            "presentation__max_categories": 100,
            "single_dataset__obs_names": None,
            "single_dataset__var_names": None,
            "diffexp__lfc_cutoff": 0.01,
        }

    def stdAsserts(self, data):
        """ run these each time we load the data """
        self.assertIsNotNone(data)
        self.assertEqual(data.cell_count, 2638)
        self.assertEqual(data.gene_count, 1838)

    def test_posix_file(self):
        locator = DataLocator("../example-dataset/pbmc3k.h5ad")
        config = AppConfig()
        config.update(**self.args)
        config.update(single_dataset__datapath=locator.path)
        config.complete_config()
        data = AnndataAdaptor(locator, config)
        self.stdAsserts(data)

    def test_url_https(self):
        url = "https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/example-dataset/pbmc3k.h5ad"
        locator = DataLocator(url)
        config = AppConfig()
        config.update(**self.args)
        data = AnndataAdaptor(locator, config)
        self.stdAsserts(data)

    def test_url_http(self):
        url = "http://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/example-dataset/pbmc3k.h5ad"
        locator = DataLocator(url)
        config = AppConfig()
        config.update(**self.args)
        data = AnndataAdaptor(locator, config)
        self.stdAsserts(data)
