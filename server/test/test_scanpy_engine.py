import json
from os import path
import pytest
import time
import unittest
import decode_fbs

import numpy as np
from pandas import Series

from server.app.scanpy_engine.scanpy_engine import ScanpyEngine
from server.app.util.errors import FilterError


class EngineTest(unittest.TestCase):
    def setUp(self):
        args = {
            "layout": "umap",
            "diffexp": "ttest",
            "max_category_items": 100,
            "obs_names": None,
            "var_names": None,
            "diffexp_lfc_cutoff": 0.01,
        }
        self.data = ScanpyEngine("example-dataset/pbmc3k.h5ad", args)

    def test_init(self):
        self.assertEqual(self.data.cell_count, 2638)
        self.assertEqual(self.data.gene_count, 1838)
        epsilon = 0.000_005
        self.assertTrue(self.data.data.X[0, 0] - -0.171_469_51 < epsilon)

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
                "var": {"index": [1, 99, [200, 300]]}
            }
        }
        fbs = self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 102)

    def test_filter_complex(self):
        filter_ = {
            "filter": {
                "var": {
                    "annotation_value": [
                        {"name": "n_cells", "min": 10}
                    ],
                    "index": [1, 99, [200, 300]]
                }
            }
        }
        fbs = self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 91)

    def test_obs_and_var_names(self):
        self.assertEqual(np.sum(self.data.data.var["name"].isna()), 0)
        self.assertEqual(np.sum(self.data.data.obs["name"].isna()), 0)

    def test_schema(self):
        with open(path.join(path.dirname(__file__), "schema.json")) as fh:
            schema = json.load(fh)
            self.assertEqual(self.data.schema, schema)

    def test_schema_produces_error(self):
        self.data.data.obs["time"] = Series(
            list([time.time() for i in range(self.data.cell_count)]),
            dtype="datetime64[ns]",
        )
        with pytest.raises(TypeError):
            self.data._create_schema()

    def test_config(self):
        self.assertEqual(
            self.data.features["layout"]["obs"],
            {"available": True, "interactiveLimit": 50000},
        )

    def test_layout(self):
        fbs = self.data.layout_to_fbs_matrix()
        layout = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(layout["n_cols"], 2)
        self.assertEqual(layout["n_rows"], 2638)

        X = layout["columns"][0]
        self.assertTrue((X >= 0).all() and (X <= 1).all())
        Y = layout["columns"][1]
        self.assertTrue((Y >= 0).all() and (Y <= 1).all())

    def test_annotations(self):
        fbs = self.data.annotation_to_fbs_matrix("obs")
        annotations = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(annotations["n_rows"], 2638)
        self.assertEqual(annotations["n_cols"], 5)
        self.assertEqual(
            annotations["col_idx"],
            ["name", "n_genes", "percent_mito", "n_counts", "louvain"],
        )

        fbs = self.data.annotation_to_fbs_matrix("var")
        annotations = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(annotations['n_rows'], 1838)
        self.assertEqual(annotations['n_cols'], 2)
        self.assertEqual(annotations["col_idx"], ["name", "n_cells"])

    def test_annotation_fields(self):
        fbs = self.data.annotation_to_fbs_matrix("obs", ["n_genes", "n_counts"])
        annotations = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(annotations["n_rows"], 2638)
        self.assertEqual(annotations['n_cols'], 2)

        fbs = self.data.annotation_to_fbs_matrix("var", ["name"])
        annotations = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(annotations['n_rows'], 1838)
        self.assertEqual(annotations['n_cols'], 1)

    def test_diffexp_topN(self):
        f1 = {"filter": {"obs": {"index": [[0, 500]]}}}
        f2 = {"filter": {"obs": {"index": [[500, 1000]]}}}
        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"]))
        self.assertEqual(len(result), 10)
        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"], 20))
        self.assertEqual(len(result), 20)

    def test_data_frame(self):
        fbs = self.data.data_frame_to_fbs_matrix(None, "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 1838)

        with self.assertRaises(ValueError):
            self.data.data_frame_to_fbs_matrix(None, "obs")

    def test_filtered_data_frame(self):
        filter_ = {
            "filter": {"var": {"annotation_value": [{"name": "n_cells", "min": 100}]}}
        }
        fbs = self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 1040)

        filter_ = {
            "filter": {"obs": {"annotation_value": [{"name": "n_counts", "min": 3000}]}}
        }
        with self.assertRaises(FilterError):
            self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")

    def test_data_named_gene(self):
        filter_ = {
            "filter": {
                "var": {"annotation_value": [{"name": "name", "values": ["RER1"]}]}
            }
        }
        fbs = self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 1)
        self.assertEqual(data["col_idx"], [4])

        filter_ = {
            "filter": {
                "var": {"annotation_value": [{"name": "name", "values": ["SPEN", "TYMP", "PRMT2"]}]}
            }
        }
        fbs = self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 3)
        self.assertTrue((data["col_idx"] == [15, 1818, 1837]).all())

    if __name__ == "__main__":
        unittest.main()
