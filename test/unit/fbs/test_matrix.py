import unittest
import pandas as pd
import numpy as np
from scipy import sparse
from parameterized import parameterized_class
import json

from test import decode_fbs
from server.common.fbs.matrix import encode_matrix_fbs, decode_matrix_fbs
from server.common.utils.type_conversion_utils import get_dtypes_and_schemas_of_dataframe
import server.common.fbs as fbs


class FbsTests(unittest.TestCase):
    """Test Case for Matrix FBS data encode/decode"""

    def test_encode_boundary(self):
        """test various boundary checks"""

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
        expected_types = ((np.ndarray, np.float32), (np.ndarray, np.int32), (np.ndarray, np.int32), (list, None))
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


"""
Test type consistency between FBS encoding and the underlying schema hint.

Basic assertion:  the FBS type returned by encode_matrix_fbs() will be consistent
with the schema hint returned by type_conversion_utils (which is in turn used
to create the client schema).

The following test cases are all dicts which contain the following keys:
    - dataframe - the dataframe used as input for encode_matrix_fbs
    - expected_fbs_types - upon success, dict of FBS column types expected (eg, Float32Array)
    - expected_schema_hints - upon success, dict of schema hint
All are keyed by column name.
"""

# simple tests that we convert all ints to int32
int_dtypes = [np.dtype(d) for d in [np.int8, np.uint8, np.int16, np.uint16, np.int32, np.uint32, np.int64, np.uint64]]
int_test_cases = [
    {
        "dataframe": pd.DataFrame({dtype.name: np.zeros((10,), dtype=dtype) for dtype in int_dtypes}),
        "expected_fbs_types": dict(
            [(dtype.name, fbs.NetEncoding.TypedArray.TypedArray.Int32Array) for dtype in int_dtypes]
        ),
        "expected_schema_hints": dict([(dtype.name, {"type": "int32"}) for dtype in int_dtypes]),
    }
]

# simple tests that we convert all floats to float32
float_dtypes = [np.dtype(d) for d in [np.float16, np.float32, np.float64]]
float_test_cases = [
    {
        "dataframe": pd.DataFrame({dtype.name: np.zeros((10,), dtype=dtype) for dtype in float_dtypes}),
        "expected_fbs_types": dict(
            [(dtype.name, fbs.NetEncoding.TypedArray.TypedArray.Float32Array) for dtype in float_dtypes]
        ),
        "expected_schema_hints": dict([(dtype.name, {"type": "float32"}) for dtype in float_dtypes]),
    }
]

# boolean - should be encoded as an uint32
bool_dtypes = [np.dtype(d) for d in [np.bool_, bool]]
bool_test_cases = [
    {
        "dataframe": pd.DataFrame({dtype.name: np.ones((10,), dtype=dtype) for dtype in bool_dtypes}),
        "expected_fbs_types": dict(
            [(dtype.name, fbs.NetEncoding.TypedArray.TypedArray.Uint32Array) for dtype in bool_dtypes]
        ),
        "expected_schema_hints": dict([(dtype.name, {"type": "boolean"}) for dtype in bool_dtypes]),
    }
]

cat_test_cases = [
    {
        "dataframe": pd.DataFrame({"a": pd.Series(["a", "b", "c", "a", "b", "c"], dtype="category")}),
        "expected_fbs_types": {"a": fbs.NetEncoding.TypedArray.TypedArray.JSONEncodedArray},
        "expected_schema_hints": {"a": {"type": "categorical", "categories": ["a", "b", "c"]}},
    },
    {
        "dataframe": pd.DataFrame(
            {"a": pd.Series(["a", "b", "c", "a", "b", "c"], dtype="category").cat.remove_categories("b")}
        ),
        "expected_fbs_types": {"a": fbs.NetEncoding.TypedArray.TypedArray.JSONEncodedArray},
        "expected_schema_hints": {"a": {"type": "categorical", "categories": ["a", "c"]}},
    },
    {
        "dataframe": pd.DataFrame({"a": pd.Series(np.arange(0, 10, dtype=np.int64), dtype="category")}),
        "expected_fbs_types": {"a": fbs.NetEncoding.TypedArray.TypedArray.Int32Array},
        "expected_schema_hints": {"a": {"type": "categorical"}},
    },
    {
        "dataframe": pd.DataFrame(
            {"a": pd.Series(np.arange(0, 10, dtype=np.int64), dtype="category").cat.remove_categories(2)}
        ),
        "expected_fbs_types": {"a": fbs.NetEncoding.TypedArray.TypedArray.Float32Array},
        "expected_schema_hints": {"a": {"type": "categorical"}},
    },
    {
        "dataframe": pd.DataFrame({"a": pd.Series(np.arange(0, 10, dtype=np.float64), dtype="category")}),
        "expected_fbs_types": {"a": fbs.NetEncoding.TypedArray.TypedArray.Float32Array},
        "expected_schema_hints": {"a": {"type": "categorical"}},
    },
    {
        "dataframe": pd.DataFrame(
            {"a": pd.Series(np.arange(0, 10, dtype=np.float64), dtype="category").cat.remove_categories(2)}
        ),
        "expected_fbs_types": {"a": fbs.NetEncoding.TypedArray.TypedArray.Float32Array},
        "expected_schema_hints": {"a": {"type": "categorical"}},
    },
]


test_cases = [
    *int_test_cases,
    *float_test_cases,
    *bool_test_cases,
    *cat_test_cases,
]


@parameterized_class(test_cases)
class TestTypeConversionConsistency(unittest.TestCase):
    def test_type_conversion_consistency(self):
        self.assertEqual(self.dataframe.shape[1], len(self.expected_fbs_types))
        self.assertEqual(self.dataframe.shape[1], len(self.expected_schema_hints))

        buf = encode_matrix_fbs(matrix=self.dataframe, col_idx=self.dataframe.columns)
        encoding_dtypes, schema_hints = get_dtypes_and_schemas_of_dataframe(self.dataframe)

        # check schema hints
        # print(schema_hints)
        # print(self.expected_schema_hints)
        self.assertEqual(schema_hints, self.expected_schema_hints)

        # inspect the FBS types
        matrix = fbs.NetEncoding.Matrix.Matrix.GetRootAsMatrix(buf, 0)
        columns_length = matrix.ColumnsLength()
        self.assertEqual(columns_length, self.dataframe.shape[1])

        self.assertEqual(matrix.ColIndexType(), fbs.NetEncoding.TypedArray.TypedArray.JSONEncodedArray)
        col_labels_arr = fbs.NetEncoding.JSONEncodedArray.JSONEncodedArray()
        col_labels_arr.Init(matrix.ColIndex().Bytes, matrix.ColIndex().Pos)
        col_index_labels = json.loads(col_labels_arr.DataAsNumpy().tobytes().decode("utf-8"))
        self.assertEqual(len(col_index_labels), self.dataframe.shape[1])

        for col_idx in range(0, columns_length):
            col_label = col_index_labels[col_idx]
            col = matrix.Columns(col_idx)
            col_type = col.UType()
            self.assertEqual(self.expected_fbs_types[col_label], col_type)
