import json
from os import path, listdir
import unittest
import server.test.decode_fbs as decode_fbs
import shutil

import numpy as np
import pandas as pd

from server.common.rest import schema_get_helper, annotations_put_fbs_helper
from server.test import data_with_tmp_annotations, make_fbs
from server.data_common.matrix_loader import MatrixDataType


class WritableAnnotationTest(unittest.TestCase):
    def setUp(self):
        self.data, self.tmp_dir, self.annotations = data_with_tmp_annotations(MatrixDataType.H5AD)
        self.data.dataset_config.user_annotations = self.annotations

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def annotation_put_fbs(self, fbs):
        annotations_put_fbs_helper(self.data, fbs)
        res = json.dumps({"status": "OK"})
        return res

    def test_error_checks(self):
        # verify that the expected errors are generated
        n_rows = self.data.get_shape()[0]
        fbs_bad = make_fbs({"louvain": pd.Series(["undefined"] * n_rows, dtype="category")})

        # ensure we catch attempt to overwrite non-writable data
        with self.assertRaises(KeyError):
            self.annotation_put_fbs(fbs_bad)

    def test_write_to_file(self):
        # verify the file is written as expected
        n_rows = self.data.get_shape()[0]
        fbs = make_fbs(
            {
                "cat_A": pd.Series(["label_A"] * n_rows, dtype="category"),
                "cat_B": pd.Series(["label_B"] * n_rows, dtype="category"),
            }
        )
        res = self.annotation_put_fbs(fbs)
        self.assertEqual(res, json.dumps({"status": "OK"}))
        self.assertTrue(path.exists(self.annotations.output_file))
        df = pd.read_csv(self.annotations.output_file, index_col=0, header=0, comment="#")
        self.assertEqual(df.shape, (n_rows, 2))
        self.assertEqual(set(df.columns), {"cat_A", "cat_B"})
        self.assertTrue(self.data.original_obs_index.equals(df.index))
        self.assertTrue(np.all(df["cat_A"] == ["label_A"] * n_rows))
        self.assertTrue(np.all(df["cat_B"] == ["label_B"] * n_rows))

        # verify complete overwrite on second attempt, AND rotation occurs
        fbs = make_fbs(
            {
                "cat_A": pd.Series(["label_A1"] * n_rows, dtype="category"),
                "cat_C": pd.Series(["label_C"] * n_rows, dtype="category"),
            }
        )
        res = self.annotation_put_fbs(fbs)
        self.assertEqual(res, json.dumps({"status": "OK"}))
        self.assertTrue(path.exists(self.annotations.output_file))
        df = pd.read_csv(self.annotations.output_file, index_col=0, header=0, comment="#")
        self.assertEqual(set(df.columns), {"cat_A", "cat_C"})
        self.assertTrue(np.all(df["cat_A"] == ["label_A1"] * n_rows))
        self.assertTrue(np.all(df["cat_C"] == ["label_C"] * n_rows))

        # rotation
        name, ext = path.splitext(self.annotations.output_file)
        backup_dir = f"{name}-backups"
        self.assertTrue(path.isdir(backup_dir))
        found_files = listdir(backup_dir)
        self.assertEqual(len(found_files), 1)

    def test_file_rotation_to_max_9(self):
        # verify we stop rotation at 9
        n_rows = self.data.get_shape()[0]
        fbs = make_fbs(
            {
                "cat_A": pd.Series(["label_A"] * n_rows, dtype="category"),
                "cat_B": pd.Series(["label_B"] * n_rows, dtype="category"),
            }
        )
        for i in range(0, 11):
            res = self.annotation_put_fbs(fbs)
            self.assertEqual(res, json.dumps({"status": "OK"}))

        name, ext = path.splitext(self.annotations.output_file)
        backup_dir = f"{name}-backups"
        self.assertTrue(path.isdir(backup_dir))
        found_files = listdir(backup_dir)
        self.assertTrue(len(found_files) <= 9)

    def test_put_get_roundtrip(self):
        # verify that OBS PUTs (annotation_put_fbs) are accessible via
        # GET (annotation_to_fbs_matrix)

        n_rows = self.data.get_shape()[0]
        fbs = make_fbs(
            {
                "cat_A": pd.Series(["label_A"] * n_rows, dtype="category"),
                "cat_B": pd.Series(["label_B"] * n_rows, dtype="category"),
            }
        )

        # put
        res = self.annotation_put_fbs(fbs)
        self.assertEqual(res, json.dumps({"status": "OK"}))

        # get
        labels = self.annotations.read_labels(None)
        fbsAll = self.data.annotation_to_fbs_matrix("obs", None, labels)
        schema = schema_get_helper(self.data)
        annotations = decode_fbs.decode_matrix_FBS(fbsAll)
        obs_index_col_name = schema["annotations"]["obs"]["index"]
        self.assertEqual(annotations["n_rows"], n_rows)
        self.assertEqual(annotations["n_cols"], 7)
        self.assertIsNone(annotations["row_idx"])
        self.assertEqual(
            annotations["col_idx"],
            [obs_index_col_name, "n_genes", "percent_mito", "n_counts", "louvain", "cat_A", "cat_B"],
        )
        col_idx = annotations["col_idx"]
        self.assertEqual(annotations["columns"][col_idx.index("cat_A")], ["label_A"] * n_rows)
        self.assertEqual(annotations["columns"][col_idx.index("cat_B")], ["label_B"] * n_rows)

        # verify the schema was updated
        all_col_schema = {c["name"]: c for c in schema["annotations"]["obs"]["columns"]}
        self.assertEqual(
            all_col_schema["cat_A"],
            {"name": "cat_A", "type": "categorical", "categories": ["label_A"], "writable": True},
        )
        self.assertEqual(
            all_col_schema["cat_B"],
            {"name": "cat_B", "type": "categorical", "categories": ["label_B"], "writable": True},
        )

    def test_config(self):
        features = self.data.get_features(self.annotations)

        # test each for singular presence and accuracy of available flag
        def check_feature(method, path, available):
            feature = list(
                filter(lambda f: f.method == method and f.path == path and f.available == available, features)
            )
            self.assertIsNotNone(feature)
            self.assertEqual(len(feature), 1)

        check_feature("POST", "/cluster/", False)
        check_feature("POST", "/diffexp/", self.data.dataset_config.diffexp__enable)
        check_feature("GET", "/layout/obs", True)
        check_feature("PUT", "/layout/obs", self.data.dataset_config.embeddings__enable_reembedding)
        check_feature("PUT", "/annotations/obs", True)
