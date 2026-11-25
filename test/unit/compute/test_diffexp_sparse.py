import unittest
from unittest.mock import patch

import numpy as np
from scipy import sparse

from server.common.compute import diffexp_generic
from server.common.constants import XApproximateDistribution as Distribution


class SparseDiffExpTest(unittest.TestCase):
    def test_matches_dense_group_summaries_and_results(self):
        rng = np.random.default_rng(124)
        dense = rng.poisson(4, (60, 21)).astype(np.float32)
        dense[rng.random(dense.shape) > 0.3] = 0
        dense[:, 1] = 3
        dense[:, 2] = 0
        for constructor in (sparse.csr_matrix, sparse.csc_matrix):
            for mode in (Distribution.COUNT, Distribution.NORMAL):
                with self.subTest(constructor=constructor, mode=mode):
                    with patch.object(diffexp_generic, "_SPARSE_BLOCK_ENTRIES", 45):
                        summaries = [
                            diffexp_generic.mean_var_n(constructor(group), mode)
                            for group in (dense[:23], dense[23:])
                        ]
                    expected = [
                        diffexp_generic.mean_var_n(group.astype(np.float64), mode)
                        for group in (dense[:23], dense[23:])
                    ]
                    for actual, reference in zip(summaries, expected):
                        self.assertEqual(actual[2], reference[2])
                        np.testing.assert_allclose(actual[0], reference[0], rtol=1e-6, atol=1e-10)
                        np.testing.assert_allclose(actual[1], reference[1], rtol=1e-6, atol=1e-10)
                    actual = diffexp_generic.diffexp_ttest_from_mean_var(*summaries[0], *summaries[1], 5, 0.01)
                    reference = diffexp_generic.diffexp_ttest_from_mean_var(*expected[0], *expected[1], 5, 0.01)
                    for direction in actual:
                        np.testing.assert_allclose(actual[direction], reference[direction], rtol=1e-6, atol=1e-10)

    def test_duplicate_entries_are_summed_before_count_transform(self):
        inputs = (
            sparse.csr_matrix(([2.0, 3.0, 0.0, 1.0], [0, 0, 1, 1], [0, 3, 4]), shape=(2, 3)),
            sparse.csc_matrix(([2.0, 3.0, 0.0, 1.0], [0, 0, 0, 1], [0, 2, 4, 4]), shape=(2, 3)),
        )
        for matrix in inputs:
            with self.subTest(format=matrix.format):
                before = matrix.data.copy()
                expected = diffexp_generic.mean_var_n(matrix.toarray(), Distribution.COUNT)
                actual = diffexp_generic.mean_var_n(matrix, Distribution.COUNT)
                np.testing.assert_array_equal(matrix.data, before)
                np.testing.assert_allclose(actual[0], expected[0], rtol=1e-12)
                np.testing.assert_allclose(actual[1], expected[1], rtol=1e-12)

    def test_nonfinite_columns_keep_full_group_size(self):
        x = np.array([[np.nan, np.inf, -1.0, 0.0], [2.0, 1.0, 2.0, 0.0], [5.0, 2.0, 4.0, 1.0]])
        for constructor in (sparse.csr_matrix, sparse.csc_matrix):
            with self.subTest(constructor=constructor):
                mean, variance, n = diffexp_generic.mean_var_n(constructor(x), Distribution.COUNT)
                self.assertEqual(n, 3)
                np.testing.assert_array_equal(mean[:3], [0, 0, 0])
                np.testing.assert_array_equal(variance[:3], [0, 0, 0])
                self.assertGreater(mean[-1], 0)
                self.assertGreater(variance[-1], 0)

    def test_equal_float64_constants_have_equal_exact_means_across_batches(self):
        for constructor in (sparse.csr_matrix, sparse.csc_matrix):
            with self.subTest(constructor=constructor), patch.object(diffexp_generic, "_SPARSE_BLOCK_ENTRIES", 6):
                for n in (10, 17):
                    mean, var, count = diffexp_generic.mean_var_n(constructor(np.full((n, 2), 0.1)))
                    self.assertEqual(n, count)
                    np.testing.assert_array_equal(mean, [0.1, 0.1])
                    np.testing.assert_array_equal(var, [0, 0])

    def test_empty_and_single_row_groups(self):
        for cls in (sparse.csr_matrix, sparse.csc_matrix):
            for n in (0, 1):
                with self.subTest(format=cls.__name__, rows=n):
                    mean, var, count = diffexp_generic.mean_var_n(cls(np.full((n, 2), 0.1)))
                    self.assertEqual(count, n)
                    np.testing.assert_array_equal(mean, [0.1 * n] * 2)
                    np.testing.assert_array_equal(var, [0, 0])
