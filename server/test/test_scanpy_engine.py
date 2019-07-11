import json
from os import path, listdir
import pytest
import time
import unittest
import decode_fbs
import tempfile
import shutil

import numpy as np
import pandas as pd

from server.app.scanpy_engine.scanpy_engine import ScanpyEngine
from server.app.util.errors import FilterError, DisabledFeatureError
from server.app.util.fbs.matrix import encode_matrix_fbs


class EngineTest(unittest.TestCase):
    def setUp(self):
        args = {
            "layout": ["umap"],
            "max_category_items": 100,
            "obs_names": None,
            "var_names": None,
            "diffexp_lfc_cutoff": 0.01,
            "layout_file": None,
        }
        self.data = ScanpyEngine("example-dataset/pbmc3k.h5ad", args)

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
        self.assertEqual(np.sum(self.data.data.var[self.data.get_schema()["annotations"]["var"]["index"]].isna()), 0)
        self.assertEqual(np.sum(self.data.data.obs[self.data.get_schema()["annotations"]["obs"]["index"]].isna()), 0)

    def test_get_schema(self):
        with open(path.join(path.dirname(__file__), "schema.json")) as fh:
            schema = json.load(fh)
            self.assertEqual(self.data.get_schema(), schema)

    def test_schema_produces_error(self):
        self.data.data.obs["time"] = pd.Series(
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
        obs_index_col_name = self.data.get_schema()["annotations"]["obs"]["index"]
        self.assertEqual(
            annotations["col_idx"],
            [obs_index_col_name, "n_genes", "percent_mito", "n_counts", "louvain"],
        )

        fbs = self.data.annotation_to_fbs_matrix("var")
        annotations = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(annotations['n_rows'], 1838)
        self.assertEqual(annotations['n_cols'], 2)
        var_index_col_name = self.data.get_schema()["annotations"]["var"]["index"]
        self.assertEqual(annotations["col_idx"], [var_index_col_name, "n_cells"])

    def test_annotation_fields(self):
        fbs = self.data.annotation_to_fbs_matrix("obs", ["n_genes", "n_counts"])
        annotations = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(annotations["n_rows"], 2638)
        self.assertEqual(annotations['n_cols'], 2)

        var_index_col_name = self.data.get_schema()["annotations"]["var"]["index"]
        fbs = self.data.annotation_to_fbs_matrix("var", [var_index_col_name])
        annotations = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(annotations['n_rows'], 1838)
        self.assertEqual(annotations['n_cols'], 1)

    def test_annotation_put(self):
        with self.assertRaises(DisabledFeatureError):
            self.data.annotation_put_fbs(None, "obs")

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
        var_index_col_name = self.data.get_schema()["annotations"]["var"]["index"]
        filter_ = {
            "filter": {
                "var": {"annotation_value": [{"name": var_index_col_name, "values": ["RER1"]}]}
            }
        }
        fbs = self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 1)
        self.assertEqual(data["col_idx"], [4])

        filter_ = {
            "filter": {
                "var": {"annotation_value": [{"name": var_index_col_name, "values": ["SPEN", "TYMP", "PRMT2"]}]}
            }
        }
        fbs = self.data.data_frame_to_fbs_matrix(filter_["filter"], "var")
        data = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(data["n_rows"], 2638)
        self.assertEqual(data["n_cols"], 3)
        self.assertTrue((data["col_idx"] == [15, 1818, 1837]).all())


class WritableAnnotationTest(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp()
        self.label_file = path.join(self.tmpDir, "labels.csv")
        args = {
            "layout": ["umap"],
            "max_category_items": 100,
            "obs_names": None,
            "var_names": None,
            "diffexp_lfc_cutoff": 0.01,
            "label_file": self.label_file
        }
        self.data = ScanpyEngine("example-dataset/pbmc3k.h5ad", args)

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def make_fbs(self, data):
        df = pd.DataFrame(data)
        return encode_matrix_fbs(matrix=df, row_idx=None, col_idx=df.columns)

    def test_error_checks(self):
        # verify that the expected errors are generated

        n_rows = self.data.data.obs.shape[0]
        fbs_bad = self.make_fbs({
            'louvain': pd.Series(['undefined' for l in range(0, n_rows)], dtype='category')
        })

        # ensure attempt to change VAR annotation
        with self.assertRaises(ValueError):
            self.data.annotation_put_fbs("var", fbs_bad)

        # ensure we catch attempt to overwrite non-writable data
        with self.assertRaises(KeyError):
            self.data.annotation_put_fbs("obs", fbs_bad)

    def test_write_to_file(self):
        # verify the file is written as expected
        n_rows = self.data.data.obs.shape[0]
        fbs = self.make_fbs({
            'cat_A': pd.Series(['label_A' for l in range(0, n_rows)], dtype='category'),
            'cat_B': pd.Series(['label_B' for l in range(0, n_rows)], dtype='category')
        })
        res = self.data.annotation_put_fbs("obs", fbs)
        self.assertEqual(res, json.dumps({"status": "OK"}))
        self.assertTrue(path.exists(self.label_file))
        df = pd.read_csv(self.label_file)
        self.assertEqual(df.shape, (n_rows, 2))
        self.assertEqual(set(df.columns), set(['cat_A', 'cat_B']))
        print(df['cat_A'].values)
        self.assertTrue(np.all(df['cat_A'] == ['label_A' for l in range(0, n_rows)]))
        self.assertTrue(np.all(df['cat_B'] == ['label_B' for l in range(0, n_rows)]))

        # verify complete overwrite on second attempt, AND rotation occurs
        fbs = self.make_fbs({
            'cat_A': pd.Series(['label_A1' for l in range(0, n_rows)], dtype='category'),
            'cat_C': pd.Series(['label_C' for l in range(0, n_rows)], dtype='category')
        })
        res = self.data.annotation_put_fbs("obs", fbs)
        self.assertEqual(res, json.dumps({"status": "OK"}))
        self.assertTrue(path.exists(self.label_file))
        df = pd.read_csv(self.label_file)
        self.assertEqual(set(df.columns), set(['cat_A', 'cat_C']))
        self.assertTrue(np.all(df['cat_A'] == ['label_A1' for l in range(0, n_rows)]))
        self.assertTrue(np.all(df['cat_C'] == ['label_C' for l in range(0, n_rows)]))

        # rotation
        name, ext = path.splitext(self.label_file)
        self.assertTrue(path.exists(f"{name}-1{ext}"))

    def test_file_rotation_to_max_9(self):
        # verify we stop rotation at 9
        n_rows = self.data.data.obs.shape[0]
        fbs = self.make_fbs({
            'cat_A': pd.Series(['label_A' for l in range(0, n_rows)], dtype='category'),
            'cat_B': pd.Series(['label_B' for l in range(0, n_rows)], dtype='category')
        })
        for i in range(0, 11):
            res = self.data.annotation_put_fbs("obs", fbs)
            self.assertEqual(res, json.dumps({"status": "OK"}))

        name, ext = path.splitext(self.label_file)
        expected_files = [self.label_file] + [f"{name}-{i}{ext}" for i in range(1, 10)]
        found_files = [path.join(self.tmpDir, p) for p in listdir(self.tmpDir)]
        self.assertEqual(set(expected_files), set(found_files))

    def test_put_get_roundtrip(self):
        # verify that OBS PUTs (annotation_put_fbs) are accessible via
        # GET (annotation_to_fbs_matrix)

        n_rows = self.data.data.obs.shape[0]
        fbs = self.make_fbs({
            'cat_A': pd.Series(['label_A' for l in range(0, n_rows)], dtype='category'),
            'cat_B': pd.Series(['label_B' for l in range(0, n_rows)], dtype='category')
        })

        # put
        res = self.data.annotation_put_fbs("obs", fbs)
        self.assertEqual(res, json.dumps({"status": "OK"}))

        # get
        fbsAll = self.data.annotation_to_fbs_matrix("obs")
        schema = self.data.get_schema()
        annotations = decode_fbs.decode_matrix_FBS(fbsAll)
        obs_index_col_name = schema["annotations"]["obs"]["index"]
        self.assertEqual(annotations["n_rows"], n_rows)
        self.assertEqual(annotations["n_cols"], 7)
        self.assertIsNone(annotations["row_idx"])
        self.assertEqual(annotations["col_idx"], [
            obs_index_col_name, "n_genes", "percent_mito", "n_counts", "louvain", "cat_A", "cat_B"
        ])
        col_idx = annotations["col_idx"]
        self.assertEqual(annotations["columns"][col_idx.index('cat_A')], [
            'label_A' for l in range(0, n_rows)
        ])
        self.assertEqual(annotations["columns"][col_idx.index('cat_B')], [
            'label_B' for l in range(0, n_rows)
        ])

        # verify the schema was updated
        all_col_schema = {c["name"]: c for c in schema["annotations"]["obs"]["columns"]}
        self.assertEqual(all_col_schema["cat_A"], {
            "name": "cat_A",
            "type": "categorical",
            "categories": ["label_A"],
            "writable": True
        })
        self.assertEqual(all_col_schema["cat_B"], {
            "name": "cat_B",
            "type": "categorical",
            "categories": ["label_B"],
            "writable": True
        })
