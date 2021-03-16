import json
import unittest
from os import path, mkdir
from shutil import rmtree
from uuid import uuid4

import numpy as np
import tiledb
from pandas import Series, DataFrame

from backend.czi_hosted.common.utils.cxg_generation_utils import (
    convert_dictionary_to_cxg_group,
    convert_dataframe_to_cxg_array,
    convert_ndarray_to_cxg_dense_array,
    convert_matrix_to_cxg_array,
)
from backend.czi_hosted.test import PROJECT_ROOT


class TestCxgGenerationUtils(unittest.TestCase):
    def setUp(self):
        self.testing_cxg_temp_directory = f"{PROJECT_ROOT}/backend/czi_hosted/test/fixtures/{uuid4()}"
        mkdir(self.testing_cxg_temp_directory)

    def tearDown(self):
        if path.isdir(self.testing_cxg_temp_directory):
            rmtree(self.testing_cxg_temp_directory)

    def test__convert_dictionary_to_cxg_group__writes_successfully(self):
        random_dictionary = {"cookies": "chocolate_chip", "brownies": "chocolate", "cake": "double chocolate"}
        dictionary_name = "favorite_desserts"
        expected_array_directory = f"{self.testing_cxg_temp_directory}/{dictionary_name}"

        convert_dictionary_to_cxg_group(
            self.testing_cxg_temp_directory, random_dictionary, group_metadata_name=dictionary_name
        )

        array = tiledb.open(expected_array_directory)
        actual_stored_metadata = dict(array.meta.items())

        self.assertTrue(path.isdir(expected_array_directory))
        self.assertTrue(isinstance(array, tiledb.DenseArray))
        self.assertEqual(random_dictionary, actual_stored_metadata)

    def test__convert_dataframe_to_cxg_array__writes_successfully(self):
        random_int_category = Series(data=[3, 1, 2, 4], dtype=np.int64)
        random_bool_category = Series(data=[True, True, False, True], dtype=np.bool_)
        random_dataframe_name = f"random_dataframe_{uuid4()}"
        random_dataframe = DataFrame(data={"int_category": random_int_category, "bool_category": random_bool_category})

        convert_dataframe_to_cxg_array(
            self.testing_cxg_temp_directory, random_dataframe_name, random_dataframe, "int_category", tiledb.Ctx()
        )

        expected_array_directory = f"{self.testing_cxg_temp_directory}/{random_dataframe_name}"
        expected_array_metadata = {
            "cxg_schema": json.dumps(
                {"int_category": {"type": "int32"}, "bool_category": {"type": "boolean"}, "index": "int_category"}
            )
        }

        actual_stored_dataframe_array = tiledb.open(expected_array_directory)
        actual_stored_dataframe_metadata = dict(actual_stored_dataframe_array.meta.items())

        self.assertTrue(path.isdir(expected_array_directory))
        self.assertTrue(isinstance(actual_stored_dataframe_array, tiledb.DenseArray))
        self.assertDictEqual(expected_array_metadata, actual_stored_dataframe_metadata)
        self.assertTrue((actual_stored_dataframe_array[0:4]["int_category"] == random_int_category.to_numpy()).all())
        self.assertTrue((actual_stored_dataframe_array[0:4]["bool_category"] == random_bool_category.to_numpy()).all())

    def test__convert_ndarray_to_cxg_dense_array__writes_successfully(self):
        ndarray = np.random.rand(3, 2)
        ndarray_name = f"{self.testing_cxg_temp_directory}/awesome_ndarray_{uuid4()}"

        convert_ndarray_to_cxg_dense_array(ndarray_name, ndarray, tiledb.Ctx())

        actual_stored_array = tiledb.open(ndarray_name)

        self.assertTrue(path.isdir(ndarray_name))
        self.assertTrue(isinstance(actual_stored_array, tiledb.DenseArray))
        self.assertTrue((actual_stored_array[:, :] == ndarray).all())

    def test__convert_matrix_to_cxg_array__dense_array_writes_successfully(self):
        matrix = np.float32(np.random.rand(3, 2))
        matrix_name = f"{self.testing_cxg_temp_directory}/awesome_matrix_{uuid4()}"

        convert_matrix_to_cxg_array(matrix_name, matrix, False, tiledb.Ctx())

        actual_stored_array = tiledb.open(matrix_name)

        self.assertTrue(path.isdir(matrix_name))
        self.assertTrue(isinstance(actual_stored_array, tiledb.DenseArray))
        self.assertTrue((actual_stored_array[:, :] == matrix).all())

    def test__convert_matrix_to_cxg_array__sparse_array_only_store_nonzeros_empty_array(self):
        matrix = np.zeros([3, 2])
        matrix_name = f"{self.testing_cxg_temp_directory}/awesome_zero_matrix_{uuid4()}"

        convert_matrix_to_cxg_array(matrix_name, matrix, True, tiledb.Ctx())

        actual_stored_array = tiledb.open(matrix_name)

        self.assertTrue(path.isdir(matrix_name))
        self.assertTrue(isinstance(actual_stored_array, tiledb.SparseArray))
        self.assertTrue(actual_stored_array[:, :][""].size == 0)

    def test__convert_matrix_to_cxg_array__sparse_array_only_store_nonzeros(self):
        matrix = np.zeros([3, 3])
        matrix[0, 0] = 1
        matrix[1, 1] = 1
        matrix[2, 2] = 2
        matrix_name = f"{self.testing_cxg_temp_directory}/awesome_sparse_matrix_{uuid4()}"

        convert_matrix_to_cxg_array(matrix_name, matrix, True, tiledb.Ctx())

        actual_stored_array = tiledb.open(matrix_name)

        self.assertTrue(path.isdir(matrix_name))
        self.assertTrue(isinstance(actual_stored_array, tiledb.SparseArray))
        self.assertTrue(actual_stored_array[0, 0][""] == 1)
        self.assertTrue(actual_stored_array[1, 1][""] == 1)
        self.assertTrue(actual_stored_array[2, 2][""] == 2)
        self.assertTrue(actual_stored_array[:, :][""].size == 3)

    def test__convert_matrix_to_cxg_array__sparse_array_with_column_encoding_empty_array(self):
        matrix_name = f"{self.testing_cxg_temp_directory}/awesome_column_shift_matrix_{uuid4()}"
        matrix = np.ones((3, 2))
        # The column shift will be equal to the matrix since subtracting the column shift from the matrix will create
        # a matrix of zeros which is sparse.
        column_shift = np.ones((3, 2))

        convert_matrix_to_cxg_array(
            matrix_name, matrix, True, tiledb.Ctx(), column_shift_for_sparse_encoding=column_shift
        )

        actual_stored_array = tiledb.open(matrix_name)

        self.assertTrue(path.isdir(matrix_name))
        self.assertTrue(isinstance(actual_stored_array, tiledb.SparseArray))
        self.assertTrue(actual_stored_array[:, :][""].size == 0)

    def test__convert_matrix_to_cxg_array__sparse_array_with_column_encoding_partial_array(self):
        matrix_name = f"{self.testing_cxg_temp_directory}/awesome_column_shift_matrix_{uuid4()}"
        matrix = np.ones((2, 2))
        # Only column shift the first column of ones.
        column_shift = np.array([[1, 0], [1, 0]])

        convert_matrix_to_cxg_array(
            matrix_name, matrix, True, tiledb.Ctx(), column_shift_for_sparse_encoding=column_shift
        )

        actual_stored_array = tiledb.open(matrix_name)

        self.assertTrue(path.isdir(matrix_name))
        self.assertTrue(isinstance(actual_stored_array, tiledb.SparseArray))
        self.assertTrue(actual_stored_array[0, 1][""] == 1)
        self.assertTrue(actual_stored_array[1, 1][""] == 1)
        self.assertTrue(actual_stored_array[:, :][""].size == 2)
