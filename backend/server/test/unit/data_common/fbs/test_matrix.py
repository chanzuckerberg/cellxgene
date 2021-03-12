import unittest
import pandas as pd
import numpy as np
from scipy import sparse

from backend.server.test.unit import decode_fbs
from backend.server.data_common.fbs.matrix import encode_matrix_fbs, decode_matrix_fbs


class FbsTests(unittest.TestCase):
    """Test Case for Matrix FBS data encode/decode """

    def test_encode_boundary(self):
        """ test various boundary checks """

        # row indexing is unsupported
        with self.assertRaises(ValueError):
            encode_matrix_fbs(matrix=pd.DataFrame(), row_idx=[])

        # matrix must be 2D
        with self.assertRaises(ValueError):
            encode_matrix_fbs(matrix=np.zeros((3, 2, 1)))
        with self.assertRaises(ValueError):
            encode_matrix_fbs(matrix=np.ones((10,)))

    def fbs_checks(self, fbs, dims, expected_types, expected_column_idx):
        d = decode_fbs.decode_matrix_FBS(fbs)
        self.assertEqual(d["n_rows"], dims[0])
        self.assertEqual(d["n_cols"], dims[1])
        self.assertIsNone(d["row_idx"])
        self.assertEqual(len(d["columns"]), dims[1])
        for i in range(0, len(d["columns"])):
            self.assertEqual(len(d["columns"][i]), dims[0])
            self.assertIsInstance(d["columns"][i], expected_types[i][0])
            if expected_types[i][1] is not None:
                self.assertEqual(d["columns"][i].dtype, expected_types[i][1])
        if expected_column_idx is not None:
            self.assertSetEqual(set(expected_column_idx), set(d["col_idx"]))

    def test_encode_DataFrame(self):
        df = pd.DataFrame(
            data={
                "a": np.zeros((10,), dtype=np.float32),
                "b": np.ones((10,), dtype=np.int64),
                "c": np.array([i for i in range(0, 10)], dtype=np.uint16),
                "d": pd.Series(["x", "y", "z", "x", "y", "z", "a", "x", "y", "z"], dtype="category"),
            }
        )
        expected_types = ((np.ndarray, np.float32), (np.ndarray, np.int32), (np.ndarray, np.uint32), (list, None))
        fbs = encode_matrix_fbs(matrix=df, row_idx=None, col_idx=df.columns)
        self.fbs_checks(fbs, (10, 4), expected_types, ["a", "b", "c", "d"])

    def test_encode_ndarray(self):
        arr = np.zeros((3, 2), dtype=np.float32)
        expected_types = ((np.ndarray, np.float32), (np.ndarray, np.float32), (np.ndarray, np.float32))
        fbs = encode_matrix_fbs(matrix=arr, row_idx=None, col_idx=None)
        self.fbs_checks(fbs, (3, 2), expected_types, None)

    def test_encode_sparse(self):
        csc = sparse.csc_matrix(np.array([[0, 1, 2], [3, 0, 4]]))
        expected_types = ((np.ndarray, np.int32), (np.ndarray, np.int32), (np.ndarray, np.int32))
        fbs = encode_matrix_fbs(matrix=csc, row_idx=None, col_idx=None)
        self.fbs_checks(fbs, (2, 3), expected_types, None)

    def test_roundtrip(self):
        dfSrc = pd.DataFrame(
            data={
                "a": np.zeros((10,), dtype=np.float32),
                "b": np.ones((10,), dtype=np.int64),
                "c": np.array([i for i in range(0, 10)], dtype=np.uint16),
                "d": pd.Series(["x", "y", "z", "x", "y", "z", "a", "x", "y", "z"], dtype="category"),
            }
        )
        dfDst = decode_matrix_fbs(encode_matrix_fbs(matrix=dfSrc, col_idx=dfSrc.columns))
        self.assertEqual(dfSrc.shape, dfDst.shape)
        self.assertEqual(set(dfSrc.columns), set(dfDst.columns))
        for c in dfSrc.columns:
            self.assertTrue(c in dfDst.columns)
            if isinstance(dfSrc[c], pd.Series):
                self.assertTrue(np.all(dfSrc[c] == dfDst[c]))
            else:
                self.assertEqual(dfSrc[c], dfDst[c])
