import json
import pytest
import unittest
import warnings
import math

import decode_fbs

from server.app.scanpy_engine.scanpy_engine import ScanpyEngine
from server.app.util.errors import JSONEncodingValueError


class NaNTest(unittest.TestCase):
    def setUp(self):
        self.args = {
            "layout": "umap",
            "diffexp": "ttest",
            "max_category_items": 100,
            "obs_names": None,
            "var_names": None,
            "diffexp_lfc_cutoff": 0.01,
        }
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            self.data = ScanpyEngine("server/test/test_datasets/nan.h5ad", self.args)
            self.data._create_schema()

    def test_load(self):
        with self.assertWarns(UserWarning):
            ScanpyEngine("server/test/test_datasets/nan.h5ad", self.args)

    def test_init(self):
        self.assertEqual(self.data.cell_count, 100)
        self.assertEqual(self.data.gene_count, 100)
        epsilon = 0.000_005
        self.assertTrue(self.data.data.X[0, 0] - -0.171_469_51 < epsilon)

    def test_dataframe(self):
        data_frame_var = decode_fbs.decode_matrix_FBS(self.data.data_frame_to_fbs_matrix(None, "var"))
        self.assertIsNotNone(data_frame_var)
        self.assertEqual(data_frame_var["n_rows"], 100)
        self.assertEqual(data_frame_var["n_cols"], 100)
        self.assertTrue(math.isnan(data_frame_var["columns"][3][3]))

        with self.assertRaises(ValueError) as cm:
            decode_fbs.decode_matrix_FBS(self.data.data_frame_to_fbs_matrix(None, "obs"))
        self.assertIsNotNone(cm.exception)

        with pytest.raises(JSONEncodingValueError):
            data_frame_obs = json.loads(self.data.data_frame(None, "obs"))
        with pytest.raises(JSONEncodingValueError):
            data_frame_var = json.loads(self.data.data_frame(None, "var"))

    def test_annotation(self):
        annotations = decode_fbs.decode_matrix_FBS(self.data.annotation_to_fbs_matrix("obs"))
        self.assertEqual(
            annotations["col_idx"],
            ["name", "n_genes", "percent_mito", "n_counts", "louvain"]
        )
        self.assertEqual(annotations["n_rows"], 100)
        self.assertTrue(math.isnan(annotations["columns"][2][0]))

        annotations = decode_fbs.decode_matrix_FBS(self.data.annotation_to_fbs_matrix("var"))
        self.assertEqual(annotations["col_idx"], ["name", "n_cells", "var_with_nans"])
        self.assertEqual(annotations["n_rows"], 100)
        self.assertTrue(math.isnan(annotations["columns"][2][0]))

        with pytest.raises(JSONEncodingValueError):
            annotations = json.loads(self.data.annotation(None, "obs"))
        with pytest.raises(JSONEncodingValueError):
            annotations = json.loads(self.data.annotation(None, "var"))
