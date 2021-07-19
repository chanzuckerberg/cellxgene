import json
import sys
import time
import unittest

import numpy as np
import pandas as pd
import pytest
from parameterized import parameterized_class

import backend.test.decode_fbs as decode_fbs
from backend.common.utils.data_locator import DataLocator
from backend.common.errors import FilterError
from backend.server.data_anndata.anndata_adaptor import AnndataAdaptor
from backend.test import PROJECT_ROOT, FIXTURES_ROOT
from backend.test.test_server.unit import app_config
from backend.test.fixtures.fixtures import pbmc3k_colors

"""
Test the anndata adaptor using the pbmc3k data set.
"""


@parameterized_class(
    ("data_locator", "backed"),
    [
        (f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad", False),
        (f"{FIXTURES_ROOT}/pbmc3k-CSC-gz.h5ad", False),
        (f"{FIXTURES_ROOT}/pbmc3k-CSR-gz.h5ad", False),
        (f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad", True),
        (f"{FIXTURES_ROOT}/pbmc3k-CSC-gz.h5ad", True),
        (f"{FIXTURES_ROOT}/pbmc3k-CSR-gz.h5ad", True),
    ],
)
class AdaptorTest(unittest.TestCase):
    def setUp(self):
        config = app_config(self.data_locator, self.backed)
        self.data = AnndataAdaptor(DataLocator(self.data_locator), config)

    def test_init(self):
        self.assertEqual(self.data.cell_count, 2638)
        self.assertEqual(self.data.gene_count, 1838)
        epsilon = 0.000_005
        self.assertTrue(self.data.data.X[0, 0] - -0.171_469_51 < epsilon)

    def test_mandatory_annotations(self):
        obs_index_col_name = self.data.get_schema()["annotations"]["obs"]["index"]
        self.assertIn(obs_index_col_name, self.data.data.obs)
        self.assertEqual(list(self.data.data.obs.index), list(range(2638)))
        var_index_col_name = self.data.get_schema()["annotations"]["var"]["index"]
        self.assertIn(var_index_col_name, self.data.data.var)
        self.assertEqual(list(self.data.data.var.index), list(range(1838)))

    @pytest.mark.filterwarnings("ignore:Anndata data matrix")
    def test_data_type(self):
        # don't run the test on the more exotic data types, as they don't
        # support the astype() interface (used by this test, but not underlying app)
        if isinstance(self.data.data.X, np.ndarray):
            self.data.data.X = self.data.data.X.astype("float64")
            with self.assertWarns(UserWarning):
                self.data._validate_data_types()

    def test_filter_idx(self):
        filter_ = {"filter": {"var": {"index": [1, 99, [200, 300]]}}}
        fbs = self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 102)

    def test_filter_complex(self):
        filter_ = {
            "filter": {"var": {"annotation_value": [{"name": "n_cells", "min": 10}], "index": [1, 99, [200, 300]]}}
        }
        fbs = self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 91)

    def test_obs_and_var_names(self):
        self.assertEqual(np.sum(self.data.data.var[self.data.get_schema()["annotations"]["var"]["index"]].isna()), 0)
        self.assertEqual(np.sum(self.data.data.obs[self.data.get_schema()["annotations"]["obs"]["index"]].isna()), 0)

    def test_get_colors(self):
        self.assertEqual(self.data.get_colors(), pbmc3k_colors)

    def test_get_schema(self):
        with open(f"{FIXTURES_ROOT}/schema.json") as fh:
            schema = json.load(fh)
            self.assertDictEqual(self.data.get_schema(), schema)

    def test_schema_produces_error(self):
        self.data.data.obs["time"] = pd.Series(
            list([time.time() for i in range(self.data.cell_count)]), dtype="datetime64[ns]",
        )
        with pytest.raises(TypeError):
            self.data._create_schema()

    def test_layout(self):
        fbs = self.data.layout_to_fbs_matrix(fields=None)
        layout = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(layout["n_cols"], 6)
        self.assertEqual(layout["n_rows"], 2638)

        X = layout["columns"][0]
        self.assertTrue((X >= 0).all() and (X <= 1).all())
        Y = layout["columns"][1]
        self.assertTrue((Y >= 0).all() and (Y <= 1).all())

    def test_layout_fields(self):
        """ X_pca, X_tsne, X_umap are available """
        fbs = self.data.layout_to_fbs_matrix(["pca"])
        layout = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(layout["n_cols"], 2)
        self.assertEqual(layout["n_rows"], 2638)
        self.assertCountEqual(layout["col_idx"], ["pca_0", "pca_1"])

        fbs = self.data.layout_to_fbs_matrix(["tsne", "pca"])
        layout = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(layout["n_cols"], 4)
        self.assertEqual(layout["n_rows"], 2638)
        self.assertCountEqual(layout["col_idx"], ["tsne_0", "tsne_1", "pca_0", "pca_1"])

    def test_annotations(self):
        fbs = self.data.annotation_to_fbs_matrix("obs")
        annotations = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(annotations["n_rows"], 2638)
        self.assertEqual(annotations["n_cols"], 5)
        obs_index_col_name = self.data.get_schema()["annotations"]["obs"]["index"]
        self.assertEqual(
            annotations["col_idx"], [obs_index_col_name, "n_genes", "percent_mito", "n_counts", "louvain"],
        )

        fbs = self.data.annotation_to_fbs_matrix("var")
        annotations = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(annotations["n_rows"], 1838)
        self.assertEqual(annotations["n_cols"], 2)
        var_index_col_name = self.data.get_schema()["annotations"]["var"]["index"]
        self.assertEqual(annotations["col_idx"], [var_index_col_name, "n_cells"])

    def test_annotation_fields(self):
        fbs = self.data.annotation_to_fbs_matrix("obs", ["n_genes", "n_counts"])
        annotations = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(annotations["n_rows"], 2638)
        self.assertEqual(annotations["n_cols"], 2)

        var_index_col_name = self.data.get_schema()["annotations"]["var"]["index"]
        fbs = self.data.annotation_to_fbs_matrix("var", [var_index_col_name])
        annotations = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(annotations["n_rows"], 1838)
        self.assertEqual(annotations["n_cols"], 1)

    def test_diffexp_topN(self):
        f1 = {"filter": {"obs": {"index": [[0, 500]]}}}
        f2 = {"filter": {"obs": {"index": [[500, 1000]]}}}
        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"]))
        self.assertEqual(len(result['positive']), 10)
        self.assertEqual(len(result['negative']), 10)

        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"], 20))
        self.assertEqual(len(result['positive']), 20)
        self.assertEqual(len(result['negative']), 20)

    def test_data_frame(self):
        f1 = {"var": {"index": [[0, 10]]}}
        fbs = self.data.data_frame_to_fbs_matrix(f1, "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 10)

        with self.assertRaises(ValueError):
            self.data.data_frame_to_fbs_matrix(None, "obs")

    def test_filtered_data_frame(self):
        filter_ = {"filter": {"var": {"annotation_value": [{"name": "n_cells", "min": 100}]}}}
        fbs = self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 1040)

        filter_ = {"filter": {"obs": {"annotation_value": [{"name": "n_counts", "min": 3000}]}}}
        with self.assertRaises(FilterError):
            self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")

    def test_data_named_gene(self):
        var_index_col_name = self.data.get_schema()["annotations"]["var"]["index"]
        filter_ = {"filter": {"var": {"annotation_value": [{"name": var_index_col_name, "values": ["RER1"]}]}}}
        fbs = self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 1)
        self.assertEqual(data["col_idx"], [4])

        filter_ = {
            "filter": {"var": {"annotation_value": [{"name": var_index_col_name, "values": ["SPEN", "TYMP", "PRMT2"]}]}}
        }
        fbs = self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 3)
        self.assertTrue((data["col_idx"] == [15, 1818, 1837]).all())
