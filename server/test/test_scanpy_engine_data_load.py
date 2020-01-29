import unittest
import json

from server.data_scanpy.scanpy_engine import ScanpyEngine
from server.common.errors import DriverError
from server.common.data_locator import DataLocator


class DataLoadEngineTest(unittest.TestCase):
    """
    Test file loading, including deferred loading/update.
    """

    def setUp(self):
        self.data_file = DataLocator("../example-dataset/pbmc3k.h5ad")
        self.data = ScanpyEngine()

    def test_init(self):
        self.assertIsNone(self.data.data)

    def test_delayed_load_args(self):
        args = {
            "layout": ["tsne"],
            "max_category_items": 1000,
            "obs_names": "foo",
            "var_names": "bar",
            "diffexp_lfc_cutoff": 0.1,
            "annotations": False,
            "annotations_file": None,
            "annotations_output_dir": None,
            "backed": False,
            "diffexp_may_be_slow": False,
            "disable_diffexp": False,
            "annotations_cell_ontology_enabled": False,
            "annotations_cell_ontology_obopath": None,
            "annotations_cell_ontology_terms": None,
        }
        self.data.update(args=args)
        self.assertEqual(args, self.data.config)

    def test_requires_data(self):
        with self.assertRaises(DriverError):
            self.data._create_schema()

    def test_delayed_load_data(self):
        self.data.update(data_locator=self.data_file)
        self.data._create_schema()
        self.assertEqual(self.data.cell_count, 2638)
        self.assertEqual(self.data.gene_count, 1838)
        epsilon = 0.000_005
        self.assertTrue(self.data.data.X[0, 0] - -0.171_469_51 < epsilon)

    def test_diffexp_topN(self):
        self.data.update(data_locator=self.data_file)
        f1 = {"filter": {"obs": {"index": [[0, 500]]}}}
        f2 = {"filter": {"obs": {"index": [[500, 1000]]}}}
        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"]))
        self.assertEqual(len(result), 10)
        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"], 20))
        self.assertEqual(len(result), 20)


class DataLocatorEngineTest(unittest.TestCase):
    """
    Test various types of data locators we expect to consume
    """

    def setUp(self):
        self.args = {
            "layout": ["umap"],
            "max_category_items": 100,
            "obs_names": None,
            "var_names": None,
            "diffexp_lfc_cutoff": 0.01,
        }

    def stdAsserts(self, data):
        """ run these each time we load the data """
        self.assertIsNotNone(data)
        self.assertEqual(data.cell_count, 2638)
        self.assertEqual(data.gene_count, 1838)

    def test_posix_file(self):
        locator = DataLocator("../example-dataset/pbmc3k.h5ad")
        data = ScanpyEngine(locator, self.args)
        self.stdAsserts(data)

    def test_url_https(self):
        url = "https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/example-dataset/pbmc3k.h5ad"
        locator = DataLocator(url)
        data = ScanpyEngine(locator, self.args)
        self.stdAsserts(data)

    def test_url_http(self):
        url = "http://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/example-dataset/pbmc3k.h5ad"
        locator = DataLocator(url)
        data = ScanpyEngine(locator, self.args)
        self.stdAsserts(data)
