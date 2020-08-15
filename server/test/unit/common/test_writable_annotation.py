import json
from os import path, listdir
import unittest
from unittest.mock import MagicMock, patch

import tiledb
from flask import Flask

import server.test.unit.decode_fbs as decode_fbs
import shutil

import numpy as np
import pandas as pd

from server.common.rest import schema_get_helper, annotations_put_fbs_helper
from server.db.cellxgene_orm import CellxGeneDataset, Annotation
from server.test import data_with_tmp_annotations, make_fbs, data_with_tmp_tiledb_annotations
from server.data_common.matrix_loader import MatrixDataType
from server.common.errors import AnnotationCategoryNameError


class auth(object):
    def get_user_id():
        return "1234"


class WritableTileDBStoredAnnotationTest(unittest.TestCase):
    def setUp(self):
        self.user_id = '1234'
        self.data, self.tmp_dir, self.annotations = data_with_tmp_tiledb_annotations(MatrixDataType.H5AD)
        self.data.dataset_config.user_annotations = self.annotations
        self.db = self.annotations.db
        self.n_rows = self.data.get_shape()[0]
        self.test_dict = {
            "cat_A": pd.Series(["label_A"] * self.n_rows, dtype="category"),
            "cat_B": pd.Series(["label_B"] * self.n_rows, dtype="category"),
        }
        self.fbs = make_fbs(self.test_dict)
        self.df = pd.DataFrame(self.test_dict)
        self.app = Flask('fake_app')
        self.app.__setattr__("auth", auth)

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def annotation_put_fbs(self, fbs):
        annotations_put_fbs_helper(self.data, fbs)
        res = json.dumps({"status": "OK"})
        return res

    def test_category_name_throws_errors_for_categories_that_cant_be_converted_to_filenames(self):
        with self.app.test_request_context():
            bad_category_names = make_fbs(
                {
                    "cat_A": pd.Series(["label_A"] * self.n_rows, dtype="category"),
                    "cat/B": pd.Series(["label_B"] * self.n_rows, dtype="category"),
                }
            )
            with self.assertRaises(AnnotationCategoryNameError):
                self.annotation_put_fbs(bad_category_names)

    def test_convert_to_pandas__converts_tiledb_to_pandas_df(self):
        with self.app.test_request_context():
            self.annotations.write_labels(self.df, self.data)
            dataset_id = self.db.query([CellxGeneDataset], [CellxGeneDataset.name == self.data.get_location()])[0].id
            annotation = self.db.query_for_most_recent(
                Annotation,
                [Annotation.user_id == self.user_id, Annotation.dataset_id == str(dataset_id)]
            )
            # retrieve tiledb array
            df = tiledb.open(annotation.tiledb_uri)
            self.assertEqual(type(df), tiledb.array.SparseArray)

            # convert to pandas df
            pandas_df = self.annotations.convert_to_pandas_df(df)
            self.assertEqual(type(pandas_df), pd.DataFrame)

    def test_write_labels_creates_a_dataset_if_it_doesnt_exist(self):
        with self.app.test_request_context():

            new_name = 'new_dataset/location'
            self.data.get_location = MagicMock(return_value=new_name)
            num_datasets = len(self.db.query([CellxGeneDataset]))
            self.annotation_put_fbs(self.fbs)
            more_datasets = len(self.db.query([CellxGeneDataset]))
            self.assertGreater(more_datasets, num_datasets)

            self.assertGreater(len(self.db.query([CellxGeneDataset], [CellxGeneDataset.name == new_name])), 0)

    def test_write_labels_links_to_existing_dataset(self):
        with self.app.test_request_context():
            # add dataset to to db
            self.annotation_put_fbs(self.fbs)

            num_datasets = len(self.db.query([CellxGeneDataset]))

            # create another annotation with the same dataset
            self.annotation_put_fbs(self.fbs)

            same_num_datasets = len(self.db.query([CellxGeneDataset]))

            self.assertEqual(num_datasets, same_num_datasets)

    def test_read_labels_returns_pandas_df(self):
        with self.app.test_request_context():
            self.annotation_put_fbs(self.fbs)
            pandas_df = self.annotations.read_labels(self.data)
            self.assertEqual(type(pandas_df), pd.DataFrame)

    def test_read_labels_returns_df_matching_original(self):
        with self.app.test_request_context():
            self.annotation_put_fbs(self.fbs)
            pandas_df = self.annotations.read_labels(self.data)

            self.assertEqual(pandas_df.shape, (self.n_rows, 2))
            self.assertEqual(set(pandas_df.columns), {"cat_A", "cat_B"})
            self.assertTrue(self.data.original_obs_index.equals(pandas_df.index))
            self.assertTrue(np.all(pandas_df["cat_A"] == ["label_A"] * self.n_rows))
            self.assertTrue(np.all(pandas_df["cat_B"] == ["label_B"] * self.n_rows))

    def test_error_checks(self):
        # verify that the expected errors are generated
        with self.app.test_request_context():
            n_rows = self.data.get_shape()[0]
            fbs_bad = make_fbs({"louvain": pd.Series(["undefined"] * n_rows, dtype="category")})

            # ensure we catch attempt to overwrite non-writable data
            with self.assertRaises(KeyError):
                self.annotation_put_fbs(fbs_bad)

    @patch('server.common.annotations.hosted_tiledb.current_app')
    def test_write_labels_stores_df_as_tiledb_array(self, mock_user_id):
        mock_user_id.auth.get_user_id.return_value = '1234'
        self.annotations.write_labels(self.df, self.data)
        # get uri
        dataset_id = self.db.query([CellxGeneDataset], [CellxGeneDataset.name == self.data.get_location()])[0].id
        annotation = self.db.query_for_most_recent(
            Annotation,
            [Annotation.user_id == '1234', Annotation.dataset_id == str(dataset_id)]
        )

        df = tiledb.open(annotation.tiledb_uri)
        self.assertEqual(type(df), tiledb.array.SparseArray)


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
