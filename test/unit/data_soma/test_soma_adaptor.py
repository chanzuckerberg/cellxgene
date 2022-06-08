import json
import time
import unittest

import numpy as np
import pandas as pd
import pytest
from parameterized import parameterized_class

import test.decode_fbs as decode_fbs
from server.common.utils.data_locator import DataLocator
from server.common.errors import FilterError
from server.data_soma.soma_adaptor import SomaAdaptor
from test import PROJECT_ROOT, FIXTURES_ROOT
from test.unit import app_config
from test.fixtures.fixtures import pbmc3k_processed_colors

"""
Test the SOMA adaptor using the pbmc3k data set.
"""


@parameterized_class(
    ("data_locator", "backed", "X_approximate_distribution"),
    [
        (f"{PROJECT_ROOT}/example-dataset/tiledb-data/pbmc3k_processed", False, "auto"),
        (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSC-gz_processed", False, "auto"),
        (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSR-gz_processed", False, "auto"),
        (f"{PROJECT_ROOT}/example-dataset/tiledb-data/pbmc3k_processed", True, "auto"),
        (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSC-gz_processed", True, "auto"),
        (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSR-gz_processed", True, "auto"),
        (f"{PROJECT_ROOT}/example-dataset/tiledb-data/pbmc3k_processed", False, "normal"),
        (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSC-gz_processed", False, "normal"),
        (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSR-gz_processed", False, "normal"),
        (f"{PROJECT_ROOT}/example-dataset/tiledb-data/pbmc3k_processed", True, "normal"),
        (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSC-gz_processed", True, "normal"),
        (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSR-gz_processed", True, "normal"),
        (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k_64_processed", False, "auto"),  # 64 bit conversion tests
        # (f"{FIXTURES_ROOT}/pbmc3k_16.h5ad", False, "auto"),  # 16 bit conversion tests are not needed for SOMA since it doesn't support 16 bit data as of now
    ],
)
class SomaAdaptorTest(unittest.TestCase):
    def setUp(self):
        config = app_config(
            self.data_locator,
            self.backed,
            extra_dataset_config=dict(X_approximate_distribution=self.X_approximate_distribution),
        )
        self.data = SomaAdaptor(DataLocator(self.data_locator), config)

    def test_init(self):
        self.assertEqual(self.data.cell_count, 2638)
        self.assertEqual(self.data.gene_count, 1838)
        epsilon = 0.000_005
        self.assertTrue(
            self.data.data.X.data.df().to_numpy()[0, 0] - -0.186_725_91 < epsilon
        )  # different hardcoded check than Anndata's test

    def test_mandatory_annotations(self):
        obs_index_col_name = self.data.get_schema()["annotations"]["obs"]["index"]
        self.assertIn(obs_index_col_name, self.data.get_obs_df())
        self.assertEqual(list(self.data.get_obs_df().index), list(range(2638)))
        var_index_col_name = self.data.get_schema()["annotations"]["var"]["index"]
        self.assertIn(var_index_col_name, self.data.get_var_df())
        self.assertEqual(list(self.data.get_var_df().index), list(range(1838)))

    # Test not needed because SOMA handles data types
    # @pytest.mark.filterwarnings("ignore:SOMA data matrix")
    # def test_data_type(self):
    #     # don't run the test on the more exotic data types, as they don't
    #     # support the astype() interface (used by this test, but not underlying app)
    #     if isinstance(self.data.data.X, np.ndarray):
    #         self.data.data.X = self.data.data.X.astype("float64")
    #         with self.assertWarns(UserWarning):
    #             self.data._validate_data_types()

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
        self.assertEqual(
            data["n_cols"], 95
        )  # different order of indices than Anndata leads to different results when using integer indices

    def test_obs_and_var_names(self):
        self.assertEqual(
            np.sum(self.data.get_var_df()[self.data.get_schema()["annotations"]["var"]["index"]].isna()), 0
        )
        self.assertEqual(
            np.sum(self.data.get_obs_df()[self.data.get_schema()["annotations"]["obs"]["index"]].isna()), 0
        )

    def test_get_colors(self):
        self.assertEqual(self.data.get_colors(), pbmc3k_processed_colors)  # different order of colors than Anndata

    def test_get_schema(self):
        with open(
            f"{FIXTURES_ROOT}/SOMA_schema.json"
        ) as fh:  # have to use a different schema than Anndata because tiledbsc return df with different datatypes
            schema = json.load(fh)
            self.assertDictEqual(self.data.get_schema(), schema)

    # This test doesn't work since TileDB-singlecell doesn't let you manually modify part of the dataset
    # def test_schema_produces_error(self):
    #     self.data.data.obs["time"] = pd.Series(
    #         list([time.time() for i in range(self.data.cell_count)]),
    #         dtype="datetime64[ns]",
    #     )
    #     with pytest.raises(TypeError):
    #         self.data._create_schema()

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
        """X_pca, X_tsne, X_umap are available"""
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
            annotations["col_idx"],
            [obs_index_col_name, "n_genes", "percent_mito", "n_counts", "louvain"],
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
        self.assertEqual(len(result["positive"]), 10)
        self.assertEqual(len(result["negative"]), 10)

        result = json.loads(self.data.diffexp_topN(f1["filter"], f2["filter"], 20))
        self.assertEqual(len(result["positive"]), 20)
        self.assertEqual(len(result["negative"]), 20)

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
        self.assertEqual(data["col_idx"], [1260])  # different than Anndata because order is not preserved

        filter_ = {
            "filter": {"var": {"annotation_value": [{"name": var_index_col_name, "values": ["SPEN", "TYMP", "PRMT2"]}]}}
        }
        fbs = self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 3)
        self.assertTrue(
            (data["col_idx"] == [1168, 1500, 1700]).all()
        )  # different than Anndata because order is not preserved
