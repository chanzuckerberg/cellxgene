import json
import pytest
import unittest
import warnings

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
            "nan_to_num": False,
        }
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            self.data = ScanpyEngine("server/test/test_datasets/nan.h5ad", self.args)
            self.data._create_schema()
        self.args_nan = dict(self.args)
        self.args_nan["nan_to_num"] = True
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            self.data_nan = ScanpyEngine(
                "server/test/test_datasets/nan.h5ad", self.args_nan
            )
            self.data_nan._create_schema()

    def test_load(self):
        with self.assertWarns(UserWarning):
            ScanpyEngine("server/test/test_datasets/nan.h5ad", self.args_nan)

    def test_init(self):
        self.assertEqual(self.data.cell_count, 100)
        self.assertEqual(self.data.gene_count, 100)
        epsilon = 0.000_005
        self.assertTrue(self.data.data.X[0, 0] - -0.171_469_51 < epsilon)

        self.assertEqual(self.data_nan.cell_count, 100)
        self.assertEqual(self.data_nan.gene_count, 100)
        epsilon = 0.000_005
        self.assertTrue(self.data_nan.data.X[0, 0] - -0.171_469_51 < epsilon)

    def test_dataframe(self):
        data_frame_obs = json.loads(self.data_nan.data_frame(None, "obs"))
        self.assertEqual(len(data_frame_obs["var"]), 100)
        self.assertEqual(len(data_frame_obs["obs"]), 100)
        data_frame_var = json.loads(self.data_nan.data_frame(None, "var"))
        self.assertEqual(len(data_frame_var["var"]), 100)
        self.assertEqual(len(data_frame_var["obs"]), 100)
        with pytest.raises(JSONEncodingValueError):
            data_frame_obs = json.loads(self.data.data_frame(None, "obs"))
        with pytest.raises(JSONEncodingValueError):
            data_frame_var = json.loads(self.data.data_frame(None, "var"))

    def test_dataframe_nan_to_0(self):
        data_frame_obs = json.loads(self.data_nan.data_frame(None, "obs"))
        self.assertEqual(data_frame_obs["obs"][1][3], 0.0)
        data_frame_var = json.loads(self.data_nan.data_frame(None, "var"))
        self.assertEqual(data_frame_var["var"][1][5], 0.0)

    def test_annotation_nan_to_0(self):
        annotations_obs = json.loads(self.data_nan.annotation(None, "obs"))
        self.assertEqual(annotations_obs["data"][0][3], 0.0)
        annotations_var = json.loads(self.data_nan.annotation(None, "var"))
        self.assertEqual(annotations_var["data"][0][3], 0.0)

    def test_annotation(self):
        annotations = json.loads(self.data_nan.annotation(None, "obs"))
        self.assertEqual(
            annotations["names"],
            ["name", "n_genes", "percent_mito", "n_counts", "louvain"],
        )
        annotations = json.loads(self.data_nan.annotation(None, "var"))
        self.assertEqual(annotations["names"], ["name", "n_cells", "var_with_nans"])
        self.assertEqual(len(annotations["data"]), 100)
        with pytest.raises(JSONEncodingValueError):
            annotations = json.loads(self.data.annotation(None, "obs"))
        with pytest.raises(JSONEncodingValueError):
            annotations = json.loads(self.data.annotation(None, "var"))
