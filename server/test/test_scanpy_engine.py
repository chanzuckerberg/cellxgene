import unittest
import pytest

from server.app.scanpy_engine.scanpy_engine import ScanpyEngine


class UtilTest(unittest.TestCase):
    def setUp(self):
        self.data = ScanpyEngine("example-dataset/")

    def test_init(self):
        self.assertEqual(self.data.cell_count, 2638)
        self.assertEqual(self.data.gene_count, 1838)
        epsilon = 0.000005
        self.assertTrue(self.data.data.X[0, 0] - -0.17146951 < epsilon)

    def test_mandatory_annotations(self):
        self.assertIn("name", self.data.data.obs)
        self.assertEqual(list(self.data.data.obs.index), list(range(2638)))
        self.assertIn("name", self.data.data.var)
        self.assertEqual(list(self.data.data.var.index), list(range(1838)))

    @pytest.mark.filterwarnings("ignore:Scanpy data matrix")
    def test_data_type(self):
        self.data.data.X = self.data.data.X.astype("float64")
        self.assertWarns(UserWarning, self.data._validatate_data_types())

    def test_filter_idx(self):
        filter_ = {
            "filter": {
                "var": {
                    "index": [1, 99, [200, 300]]
                },
                "obs": {
                    "index": [1, 99, [1000, 2000]]
                }
            }
        }
        data = self.data.filter_dataframe(filter_["filter"])
        self.assertEqual(data.shape, (1002, 102))

    def test_filter_annotation(self):
        filter_ = {
            "filter": {
                "obs": {
                    "annotation_value": [
                        {"name": "louvain", "values": ["NK cells", "CD8 T cells"]},
                    ]
                }
            }
        }
        data = self.data.filter_dataframe(filter_["filter"])
        self.assertEqual(data.shape, (470, 1838))
        filter_ = {
            "filter": {
                "obs": {
                    "annotation_value": [
                        {"name": "n_counts", "min": 3000},
                    ]
                }
            }
        }
        data = self.data.filter_dataframe(filter_["filter"])
        self.assertEqual(data.shape, (497, 1838))

    def test_filter_complex(self):
        filter_ = {
            "filter": {
                "var": {
                    "index": [1, 99, [200, 300]]
                },
                "obs": {
                    "annotation_value": [
                        {"name": "louvain", "values": ["NK cells", "CD8 T cells"]},
                        {"name": "n_counts", "min": 3000},
                    ],
                    "index": [1, 99, [1000, 2000]]
                }
            }
        }
        data = self.data.filter_dataframe(filter_["filter"])
        self.assertEqual(data.shape, (15, 102))

    if __name__ == '__main__':
        unittest.main()
