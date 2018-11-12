import json
from os import path
import pytest
import time
import unittest
import argparse

import numpy as np
from pandas import Series

from server.app.scanpy_engine.scanpy_engine import ScanpyEngine


class UtilTest(unittest.TestCase):
    def setUp(self):
        args = {'layout': 'umap', 'diffexp': 'ttest', 'max_category_items': 100,
                'obs_names': None, 'var_names': None}

        self.data = ScanpyEngine("example-dataset/pbmc3k.h5ad", args)
        self.data._create_schema()

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
        with self.assertWarns(UserWarning):
            self.data._validate_data_types()

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

    def test_filter_annotation_no_uns(self):
        filter_ = {
            "filter": {
                "var": {
                    "annotation_value": [
                        {"name": "name", "values": ["RER1"]},
                    ]
                }
            }
        }
        data = self.data.filter_dataframe(filter_["filter"], include_uns=False)
        self.assertEqual(data.shape[1], 1)

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

    def test_obs_and_var_names(self):
        self.assertEqual(np.sum(self.data.data.var["name"].isna()), 0)
        self.assertEqual(np.sum(self.data.data.obs["name"].isna()), 0)

    def test_schema(self):
        with open(path.join(path.dirname(__file__), "schema.json")) as fh:
            schema = json.load(fh)
            self.assertEqual(self.data.schema, schema)

    def test_schema_produces_error(self):
        self.data.data.obs["time"] = Series(list([time.time() for i in range(self.data.cell_count)]),
                                            dtype="datetime64[ns]")
        with pytest.raises(TypeError):
            self.data._create_schema()

    def test_config(self):
        self.assertEqual(self.data.features["layout"]["obs"], {'available': True, 'interactiveLimit': 50000})

    def test_layout(self):
        layout = self.data.layout(None)
        self.assertEqual(layout["ndims"], 2)
        self.assertEqual(len(layout["coordinates"]), 2638)
        self.assertEqual(layout["coordinates"][0][0], 0)
        for idx, val in enumerate(layout["coordinates"]):
            self.assertLessEqual(val[1], 1)
            self.assertLessEqual(val[2], 1)

    def test_annotations(self):
        annotations = self.data.annotation(None, "obs")
        self.assertEqual(annotations["names"], ["name", "n_genes", "percent_mito", "n_counts", "louvain"])
        self.assertEqual(len(annotations["data"]), 2638)
        annotations = self.data.annotation(None, "var")
        self.assertEqual(annotations["names"], ["name", "n_cells"])
        self.assertEqual(len(annotations["data"]), 1838)

    def test_annotation_fields(self):
        annotations = self.data.annotation(None, "obs", ["n_genes", "n_counts"])
        self.assertEqual(annotations["names"], ["n_genes", "n_counts"])
        self.assertEqual(len(annotations["data"]), 2638)
        annotations = self.data.annotation(None, "var", ["name"])
        self.assertEqual(annotations["names"], ["name"])
        self.assertEqual(len(annotations["data"]), 1838)

    def test_filtered_annotation(self):
        filter_ = {
            "filter": {
                "obs": {
                    "annotation_value": [
                        {"name": "n_counts", "min": 3000},
                    ]
                },
                "var": {
                    "annotation_value": [
                        {"name": "name", "values": ["ATAD3C", "RER1"]},
                    ]
                }
            }
        }
        annotations = self.data.annotation(filter_["filter"], "obs")
        self.assertEqual(annotations["names"], ["name", "n_genes", "percent_mito", "n_counts", "louvain"])
        self.assertEqual(len(annotations["data"]), 497)
        annotations = self.data.annotation(filter_["filter"], "var")
        self.assertEqual(annotations["names"], ["name", "n_cells"])
        self.assertEqual(len(annotations["data"]), 2)

    def test_filtered_layout(self):
        filter_ = {
            "filter": {
                "obs": {
                    "annotation_value": [
                        {"name": "n_counts", "min": 3000},
                    ]
                }
            }
        }
        layout = self.data.layout(filter_["filter"])
        self.assertEqual(len(layout["coordinates"]), 497)

    def test_diffexp(self):
        f1 = {
            "filter": {
                "obs": {
                    "index": [[0, 500]]
                }
            }
        }
        f2 = {
            "filter": {
                "obs": {
                    "index": [[500, 1000]]
                }
            }
        }
        result = self.data.diffexp(f1["filter"], f2["filter"])
        self.assertEqual(len(result), 10)
        var_idx = [i[0] for i in result]
        self.assertEqual(var_idx, sorted(var_idx))
        result = self.data.diffexp(f1["filter"], f2["filter"], 20)
        self.assertEqual(len(result), 20)

    def test_data_frame(self):
        data_frame_obs = self.data.data_frame(None, "obs")
        self.assertEqual(len(data_frame_obs["var"]), 1838)
        self.assertEqual(len(data_frame_obs["obs"]), 2638)
        data_frame_var = self.data.data_frame(None, "var")
        self.assertEqual(len(data_frame_var["var"]), 1838)
        self.assertEqual(len(data_frame_var["obs"]), 2638)

    def test_filtered_data_frame(self):
        filter_ = {
            "filter": {
                "obs": {
                    "annotation_value": [
                        {"name": "n_counts", "min": 3000},
                    ]
                }
            }
        }
        data_frame_obs = self.data.data_frame(filter_["filter"], "obs")
        self.assertEqual(len(data_frame_obs["var"]), 1838)
        self.assertEqual(len(data_frame_obs["obs"]), 497)
        self.assertIsInstance(data_frame_obs["obs"][0], (list, tuple))
        self.assertEqual(type(data_frame_obs["var"][0]), int)
        data_frame_var = self.data.data_frame(filter_["filter"], "var")
        self.assertEqual(len(data_frame_var["var"]), 1838)
        self.assertEqual(len(data_frame_var["obs"]), 497)
        self.assertIsInstance(data_frame_var["var"][0], (list, tuple))
        self.assertEqual(type(data_frame_var["obs"][0]), int)

    def test_data_single_gene(self):
        for axis in ["obs", "var"]:
            filter_ = {
                "filter": {
                    "var": {
                        "annotation_value": [
                            {"name": "name", "values": ["RER1"]},
                        ]
                    }
                }
            }
            data_frame_var = self.data.data_frame(filter_["filter"], axis)
            if axis == "obs":
                self.assertEqual(type(data_frame_var["var"][0]), int)
                self.assertIsInstance(data_frame_var["obs"][0], (list, tuple))
            elif axis == "var":
                self.assertEqual(type(data_frame_var["obs"][0]), int)
                self.assertIsInstance(data_frame_var["var"][0], (list, tuple))

    if __name__ == '__main__':
        unittest.main()
