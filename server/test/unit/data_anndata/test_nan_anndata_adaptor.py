import math
import unittest
import warnings

import pytest

import server.test.unit.decode_fbs as decode_fbs
from server.common.data_locator import DataLocator
from server.common.errors import FilterError
from server.data_anndata.anndata_adaptor import AnndataAdaptor
from server.test import app_config, FIXTURES_ROOT


class NaNTest(unittest.TestCase):
    def setUp(self):
        self.data_locator = DataLocator(f"{FIXTURES_ROOT}/nan.h5ad")
        self.config = app_config(self.data_locator.path)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            self.data = AnndataAdaptor(self.data_locator, self.config)
            self.data._create_schema()

    def test_load(self):
        with self.assertLogs(level="WARN") as logger:
            self.data = AnndataAdaptor(self.data_locator, self.config)
            self.assertTrue(logger.output)

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

        with pytest.raises(FilterError):
            self.data.data_frame_to_fbs_matrix("an erroneous filter", "var")
        with pytest.raises(FilterError):
            filter_ = {"filter": {"obs": {"index": [1, 99, [200, 300]]}}}
            self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")

    def test_dataframe_obs_not_implemented(self):
        with self.assertRaises(ValueError) as cm:
            decode_fbs.decode_matrix_FBS(self.data.data_frame_to_fbs_matrix(None, "obs"))
        self.assertIsNotNone(cm.exception)

    def test_annotation(self):
        annotations = decode_fbs.decode_matrix_FBS(self.data.annotation_to_fbs_matrix("obs"))
        obs_index_col_name = self.data.schema["annotations"]["obs"]["index"]
        self.assertEqual(annotations["col_idx"], [obs_index_col_name, "n_genes", "percent_mito", "n_counts", "louvain"])
        self.assertEqual(annotations["n_rows"], 100)
        self.assertTrue(math.isnan(annotations["columns"][2][0]))

        annotations = decode_fbs.decode_matrix_FBS(self.data.annotation_to_fbs_matrix("var"))
        var_index_col_name = self.data.schema["annotations"]["var"]["index"]
        self.assertEqual(annotations["col_idx"], [var_index_col_name, "n_cells", "var_with_nans"])
        self.assertEqual(annotations["n_rows"], 100)
        self.assertTrue(math.isnan(annotations["columns"][2][0]))
