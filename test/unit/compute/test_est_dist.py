import unittest
import numpy as np
from scipy import sparse
from server.common.compute.estimate_distribution import estimate_approximate_distribution
from server.common.constants import XApproximateDistribution
from server.data_common.matrix_loader import MatrixDataLoader
from test.unit import app_config
from test import PROJECT_ROOT


class EstDistTest(unittest.TestCase):
    """Tests the diffexp returns the expected results for one test case, using the h5ad
    adaptor types and different algorithms."""

    def load_dataset(self, path, extra_server_config={}, extra_dataset_config={}):
        config = app_config(path, extra_server_config=extra_server_config, extra_dataset_config=extra_dataset_config)
        loader = MatrixDataLoader(path)
        adaptor = loader.open(config)
        return adaptor

    def test_adaptestimate_approximate_distribution(self):
        adaptor = self.load_dataset(f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad")
        self.assertEqual(adaptor.get_X_approximate_distribution(), XApproximateDistribution.NORMAL)

    def test_estimate_approximate_distribution(self):
        raw = np.random.exponential(scale=1000, size=(100, 40))

        # empty
        self.assertEqual(estimate_approximate_distribution(np.zeros((0,))), XApproximateDistribution.NORMAL)

        # ndarray
        self.assertEqual(estimate_approximate_distribution(raw), XApproximateDistribution.COUNT)
        self.assertEqual(estimate_approximate_distribution(np.log1p(raw)), XApproximateDistribution.NORMAL)

        # csr_matrix
        self.assertEqual(estimate_approximate_distribution(sparse.csr_matrix(raw)), XApproximateDistribution.COUNT)
        self.assertEqual(
            estimate_approximate_distribution(sparse.csr_matrix(np.log1p(raw))), XApproximateDistribution.NORMAL
        )

        # csc_matrix
        self.assertEqual(estimate_approximate_distribution(sparse.csc_matrix(raw)), XApproximateDistribution.COUNT)
        self.assertEqual(
            estimate_approximate_distribution(sparse.csc_matrix(np.log1p(raw))), XApproximateDistribution.NORMAL
        )

        # BIG (ie, trigger MT)
        big = np.random.exponential(scale=100, size=(1_000_000, 100))
        self.assertEqual(estimate_approximate_distribution(big), XApproximateDistribution.COUNT)
        self.assertEqual(estimate_approximate_distribution(np.log1p(big)), XApproximateDistribution.NORMAL)

    def test_unsupported_throws(self):
        # dtypes and matrix formats we do not support
        with self.assertRaises(TypeError):
            estimate_approximate_distribution(np.array(["a", "b"]))
        with self.assertRaises(TypeError):
            estimate_approximate_distribution(sparse.coo_matrix(np.array([[0, 1, 2], [3, 0, 2]])))

    def test_nonfinites(self):
        def put(arr, ind, vals):
            # like np.put, but creates and returns a modified copy of original array
            a = arr.copy()
            np.put(a, ind, vals)
            return a

        # non-finites
        self.assertEqual(estimate_approximate_distribution(np.array([np.nan])), XApproximateDistribution.NORMAL)
        self.assertEqual(estimate_approximate_distribution(np.array([np.PINF])), XApproximateDistribution.NORMAL)
        self.assertEqual(estimate_approximate_distribution(np.array([np.NINF])), XApproximateDistribution.NORMAL)
        self.assertEqual(
            estimate_approximate_distribution(np.array([np.PINF, np.NINF, 0])), XApproximateDistribution.NORMAL
        )
        self.assertEqual(
            estimate_approximate_distribution(np.array([np.nan, np.PINF, np.NINF])), XApproximateDistribution.NORMAL
        )

        raw = np.random.exponential(scale=1000, size=(50, 3))
        logged = np.log1p(raw)

        self.assertEqual(
            estimate_approximate_distribution(put(raw, [1], [np.nan])),
            XApproximateDistribution.COUNT,
        )
        self.assertEqual(
            estimate_approximate_distribution(put(raw, [1], [np.PINF])),
            XApproximateDistribution.COUNT,
        )
        self.assertEqual(
            estimate_approximate_distribution(put(raw, [1], [np.NINF])),
            XApproximateDistribution.COUNT,
        )
        self.assertEqual(
            estimate_approximate_distribution(put(raw, [1, 3, 88], [np.nan, np.PINF, np.NINF])),
            XApproximateDistribution.COUNT,
        )
        self.assertEqual(
            estimate_approximate_distribution(put(raw, [0, 1], [np.nan, np.nan])),
            XApproximateDistribution.COUNT,
        )

        self.assertEqual(
            estimate_approximate_distribution(put(logged, [1], [np.nan])),
            XApproximateDistribution.NORMAL,
        )
        self.assertEqual(
            estimate_approximate_distribution(put(logged, [1], [np.PINF])),
            XApproximateDistribution.NORMAL,
        )
        self.assertEqual(
            estimate_approximate_distribution(put(logged, [1], [np.NINF])),
            XApproximateDistribution.NORMAL,
        )
        self.assertEqual(
            estimate_approximate_distribution(put(logged, [1, 3, 88], [np.nan, np.PINF, np.NINF])),
            XApproximateDistribution.NORMAL,
        )
        self.assertEqual(
            estimate_approximate_distribution(put(logged, [0, 1], [np.nan, np.nan])),
            XApproximateDistribution.NORMAL,
        )
