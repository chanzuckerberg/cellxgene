import unittest

import numpy as np

from backend.czi_hosted.common.utils.matrix_utils import is_matrix_sparse, get_column_shift_encode_for_matrix


class TestMatrixUtils(unittest.TestCase):
    def test__is_matrix_sparse__zero_and_one_hundred_percent_threshold(self):
        matrix = np.array([1, 2, 3])

        self.assertFalse(is_matrix_sparse(matrix, 0))
        self.assertTrue(is_matrix_sparse(matrix, 100))

    def test__is_matrix_sparse__partially_populated_sparse_matrix_returns_true(self):
        matrix = np.zeros([3, 4])
        matrix[2][3] = 1.0
        matrix[1][1] = 2.2

        self.assertTrue(is_matrix_sparse(matrix, 50))

    def test__is_matrix_sparse__partially_populated_dense_matrix_returns_false(self):
        matrix = np.zeros([2, 2])
        matrix[0][0] = 1.0
        matrix[0][1] = 2.2
        matrix[1][1] = 3.7

        self.assertFalse(is_matrix_sparse(matrix, 50))

    def test__is_matrix_sparse__giant_matrix_returns_false_early(self):
        matrix = np.ones([20000, 20])

        with self.assertLogs(level="INFO") as logger:
            self.assertFalse(is_matrix_sparse(matrix, 1))

            # Because the function returns early a log will output the _estimate_ instead of the _exact_ percentage of
            # non-zero elements in the matrix.
            self.assertIn("Percentage of non-zero elements (estimate)", logger.output[0])

    def test__is_matrix_sparse_with_column_shift_encoding__regular_sparse_returns_true(self):
        matrix = np.zeros([2, 2])
        matrix[0][0] = 1.0

        self.assertIsNotNone(get_column_shift_encode_for_matrix(matrix, 50))

    def test__is_matrix_sparse_with_column_shift_encoding__column_shift_returns_same_value(self):
        matrix = np.ones([2, 2])
        expected_column_shift = [1, 1]

        actual_column_shift = get_column_shift_encode_for_matrix(matrix, 50)
        self.assertTrue((expected_column_shift == actual_column_shift).all())

    def test__is_matrix_sparse_with_column_shift_encoding__impossible_column_shift_returns_none(self):
        matrix = np.array([[1, 2], [3, 4]])

        self.assertIsNone(get_column_shift_encode_for_matrix(matrix, 50))

    def test__is_matrix_sparse_with_column_shift_encoding__giant_matrix_returns_false_early(self):
        matrix = np.random.rand(20000, 20)

        with self.assertLogs(level="INFO") as logger:
            self.assertFalse(is_matrix_sparse(matrix, 1))

            # Because the function returns early a log will output the _estimate_ instead of the _exact_ percentage of
            # non-zero elements in the matrix.
            self.assertIn("Percentage of non-zero elements (estimate)", logger.output[0])
