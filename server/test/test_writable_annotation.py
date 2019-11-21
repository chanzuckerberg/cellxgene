import json
from os import path, listdir
import unittest
import decode_fbs
import tempfile
import shutil

import numpy as np
import pandas as pd

from server.app.scanpy_engine.scanpy_engine import ScanpyEngine
from server.app.util.fbs.matrix import encode_matrix_fbs
from server.app.util.data_locator import DataLocator


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
        self.data = ScanpyEngine(DataLocator("example-dataset/pbmc3k.h5ad"), args)

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
        df = pd.read_csv(self.label_file, index_col=0, header=0, comment='#')
        self.assertEqual(df.shape, (n_rows, 2))
        self.assertEqual(set(df.columns), set(['cat_A', 'cat_B']))
        self.assertTrue(self.data.original_obs_index.equals(df.index))
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
        df = pd.read_csv(self.label_file, index_col=0, header=0, comment='#')
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
