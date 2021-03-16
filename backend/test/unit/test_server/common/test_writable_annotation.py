import json
import shutil
import unittest
from os import path, listdir

import numpy as np
import pandas as pd

import backend.test.unit.test_server.decode_fbs as decode_fbs
from backend.server.common.rest import annotations_put_fbs_helper, schema_get_helper
from backend.server.data_common.matrix_loader import MatrixDataType
from backend.test.unit.test_server import data_with_tmp_annotations, make_fbs


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
        self.assertTrue(path.exists(self.annotations.label_output_file))
        df = pd.read_csv(self.annotations.label_output_file, index_col=0, header=0, comment="#")
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
        self.assertTrue(path.exists(self.annotations.label_output_file))
        df = pd.read_csv(self.annotations.label_output_file, index_col=0, header=0, comment="#")
        self.assertEqual(set(df.columns), {"cat_A", "cat_C"})
        self.assertTrue(np.all(df["cat_A"] == ["label_A1"] * n_rows))
        self.assertTrue(np.all(df["cat_C"] == ["label_C"] * n_rows))

        # rotation
        name, ext = path.splitext(self.annotations.label_output_file)
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

        name, ext = path.splitext(self.annotations.label_output_file)
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

    def test_put_float_data(self):
        # verify that OBS PUTs (annotation_put_fbs) are accessible via
        # GET (annotation_to_fbs_matrix)

        n_rows = self.data.get_shape()[0]

        # verifies that floating point with decimals fail.
        fbs = make_fbs({"cat_F_FAIL": pd.Series([1.1] * n_rows, dtype=np.dtype("float"))})
        with self.assertRaises(ValueError) as exception_context:
            res = self.annotation_put_fbs(fbs)
        self.assertEqual(str(exception_context.exception), "Columns may not have floating point types")

        # verifies that floating point that can be converted to int passes
        fbs = make_fbs({"cat_F_PASS": pd.Series([1.0] * n_rows, dtype="float")})
        res = self.annotation_put_fbs(fbs)
        self.assertEqual(res, json.dumps({"status": "OK"}))

        # check read_labels
        labels = self.annotations.read_labels(None)
        fbsAll = self.data.annotation_to_fbs_matrix("obs", None, labels)
        schema = schema_get_helper(self.data)
        annotations = decode_fbs.decode_matrix_FBS(fbsAll)
        self.assertEqual(annotations["n_rows"], n_rows)
        all_col_schema = {c["name"]: c for c in schema["annotations"]["obs"]["columns"]}
        self.assertEqual(
            all_col_schema["cat_F_PASS"],
            {"name": "cat_F_PASS", "type": "int32", "writable": True},
        )
