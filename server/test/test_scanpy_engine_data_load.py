import unittest
import json

from server.app.scanpy_engine.scanpy_engine import ScanpyEngine
from server.app.util.errors import DriverError
from server.app.util.data_locator import DataLocator


class DataLoadEngineTest(unittest.TestCase):
    def setUp(self):
        self.data_file = DataLocator("example-dataset/pbmc3k.h5ad")
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
        }
        self.data.update(args=args)
        self.assertEqual(args, self.data.config)

    def test_requires_data(self):
        with self.assertRaises(DriverError):
            self.data._create_schema()

    def test_delayed_load_data(self):
        self.data.update(data=self.data_file)
        self.data._create_schema()
        self.assertEqual(self.data.cell_count, 2638)
        self.assertEqual(self.data.gene_count, 1838)
        epsilon = 0.000_005
        self.assertTrue(self.data.data.X[0, 0] - -0.171_469_51 < epsilon)

    def test_diffexp_topN(self):
        self.data.update(data=self.data_file)
        f1 = {"filter": {"obs": {"index": [[0, 500]]}}}
        f2 = {"filter": {"obs": {"index": [[500, 1000]]}}}
        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"]))
        self.assertEqual(len(result), 10)
        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"], 20))
        self.assertEqual(len(result), 20)

    if __name__ == "__main__":
        unittest.main()
