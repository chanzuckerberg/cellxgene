import json

import numpy as np
import pandas as pd
from flatbuffers import Builder
from scipy import sparse

from server.common.utils.type_conversion_utils import get_encoding_dtype_of_array

import server.common.fbs.NetEncoding.Column as Column
import server.common.fbs.NetEncoding.Float32Array as Float32Array
import server.common.fbs.NetEncoding.Float64Array as Float64Array
import server.common.fbs.NetEncoding.Int32Array as Int32Array
import server.common.fbs.NetEncoding.JSONEncodedArray as JSONEncodedArray
import server.common.fbs.NetEncoding.Matrix as Matrix
import server.common.fbs.NetEncoding.TypedArray as TypedArray
import server.common.fbs.NetEncoding.Uint32Array as Uint32Array


# Serialization helper
def serialize_column(builder, typed_arr):
    """ Serialize NetEncoding.Column """

    (u_type, u_value) = typed_arr
    Column.ColumnStart(builder)
    Column.ColumnAddUType(builder, u_type)
    Column.ColumnAddU(builder, u_value)
    return Column.ColumnEnd(builder)


# Serialization helper
def serialize_matrix(builder, n_rows, n_cols, columns, col_idx):
    """ Serialize NetEncoding.Matrix """

    Matrix.MatrixStart(builder)
    Matrix.MatrixAddNRows(builder, n_rows)
    Matrix.MatrixAddNCols(builder, n_cols)
    Matrix.MatrixAddColumns(builder, columns)
    if col_idx is not None:
        (u_type, u_val) = col_idx
        Matrix.MatrixAddColIndexType(builder, u_type)
        Matrix.MatrixAddColIndex(builder, u_val)
    return Matrix.MatrixEnd(builder)


# Serialization helper
def serialize_typed_array(builder, source_array, encoding_info):
    """
    Serialize any of the various typed arrays, eg, Float32Array. Specific means of serialization and type conversion
    are provided by type_info.
    """

    arr = source_array
    (array_type, as_type) = encoding_info(source_array)

    if isinstance(arr, pd.Index):
        arr = arr.to_series()

    # convert to a simple ndarray
    if as_type == "json":
        as_json = arr.to_json(orient="records")
        arr = np.array(bytearray(as_json, "utf-8"))
    else:
        if sparse.issparse(arr):
            arr = arr.toarray()
        elif isinstance(arr, pd.Series):
            arr = arr.to_numpy()
        if arr.dtype != as_type:
            arr = arr.astype(as_type)

    # serialize the ndarray into a vector
    if arr.ndim == 2:
        if arr.shape[0] == 1:
            arr = arr[0]
        elif arr.shape[1] == 1:
            arr = arr.T[0]

    vec = builder.CreateNumpyVector(arr)

    # serialize the typed array table
    builder.StartObject(1)
    builder.PrependUOffsetTRelativeSlot(0, vec, 0)
    array_value = builder.EndObject()
    return (array_type, array_value)


def column_encoding(arr):
    column_encoding_type_map = {
        # array protocol string:  ( array_type, as_type )
        np.dtype(np.float64).str: (TypedArray.TypedArray.Float32Array, np.float32),
        np.dtype(np.float32).str: (TypedArray.TypedArray.Float32Array, np.float32),
        np.dtype(np.float16).str: (TypedArray.TypedArray.Float32Array, np.float32),
        np.dtype(np.int8).str: (TypedArray.TypedArray.Int32Array, np.int32),
        np.dtype(np.int16).str: (TypedArray.TypedArray.Int32Array, np.int32),
        np.dtype(np.int32).str: (TypedArray.TypedArray.Int32Array, np.int32),
        np.dtype(np.int64).str: (TypedArray.TypedArray.Int32Array, np.int32),
        np.dtype(np.uint8).str: (TypedArray.TypedArray.Uint32Array, np.uint32),
        np.dtype(np.uint16).str: (TypedArray.TypedArray.Uint32Array, np.uint32),
        np.dtype(np.uint32).str: (TypedArray.TypedArray.Uint32Array, np.uint32),
        np.dtype(np.uint64).str: (TypedArray.TypedArray.Uint32Array, np.uint32),
    }
    column_encoding_default = (TypedArray.TypedArray.JSONEncodedArray, "json")

    encoding_dtype = np.dtype(get_encoding_dtype_of_array(arr))
    return column_encoding_type_map.get(encoding_dtype.str, column_encoding_default)


def index_encoding(arr):
    index_encoding_type_map = {
        # array protocol string:  ( array_type, as_type )
        np.dtype(np.int32).str: (TypedArray.TypedArray.Int32Array, np.int32),
        np.dtype(np.int64).str: (TypedArray.TypedArray.Int32Array, np.int32),
        np.dtype(np.uint32).str: (TypedArray.TypedArray.Uint32Array, np.uint32),
        np.dtype(np.uint64).str: (TypedArray.TypedArray.Uint32Array, np.uint32),
    }
    index_encoding_default = (TypedArray.TypedArray.JSONEncodedArray, "json")

    return index_encoding_type_map.get(arr.dtype.str, index_encoding_default)


def guess_at_mem_needed(matrix):
    (n_rows, n_cols) = matrix.shape
    if isinstance(matrix, np.ndarray) or sparse.issparse(matrix):
        guess = (n_rows * n_cols * matrix.dtype.itemsize) + 1024
    elif isinstance(matrix, pd.DataFrame):
        # XXX TODO - DataFrame type estimate
        guess = 1
    else:
        guess = 1

    # round up to nearest 1024 bytes
    guess = (guess + 0x400) & (~0x3FF)
    return guess


def encode_matrix_fbs(matrix, row_idx=None, col_idx=None):
    """
    Given a 2D DataFrame, ndarray or sparse equivalent, create and return a Matrix flatbuffer.

    :param matrix: 2D DataFrame, ndarray or sparse equivalent
    :param row_idx: index for row dimension, Index or ndarray
    :param col_idx: index for col dimension, Index or ndarray

    NOTE: row indices are (currently) unsupported and must be None
    """

    if row_idx is not None:
        raise ValueError("row indexing not supported for FBS Matrix")
    if matrix.ndim != 2:
        raise ValueError("FBS Matrix must be 2D")

    (n_rows, n_cols) = matrix.shape

    # estimate size needed, so we don't unnecessarily realloc.
    builder = Builder(guess_at_mem_needed(matrix))

    columns = []
    for cidx in range(n_cols - 1, -1, -1):
        # serialize the typed array
        col = matrix.iloc[:, cidx] if isinstance(matrix, pd.DataFrame) else matrix[:, cidx]
        typed_arr = serialize_typed_array(builder, col, column_encoding)

        # serialize the Column union
        columns.append(serialize_column(builder, typed_arr))

    # Serialize Matrix.columns[]
    Matrix.MatrixStartColumnsVector(builder, n_cols)
    for c in columns:
        builder.PrependUOffsetTRelative(c)
    matrix_column_vec = builder.EndVector(n_cols)

    # serialize the colIndex if provided
    cidx = None
    if col_idx is not None:
        cidx = serialize_typed_array(builder, col_idx, index_encoding)

    # Serialize Matrix
    matrix = serialize_matrix(builder, n_rows, n_cols, matrix_column_vec, cidx)

    builder.Finish(matrix)
    return builder.Output()


def deserialize_typed_array(tarr):
    type_map = {
        TypedArray.TypedArray.NONE: None,
        TypedArray.TypedArray.Uint32Array: Uint32Array.Uint32Array,
        TypedArray.TypedArray.Int32Array: Int32Array.Int32Array,
        TypedArray.TypedArray.Float32Array: Float32Array.Float32Array,
        TypedArray.TypedArray.Float64Array: Float64Array.Float64Array,
        TypedArray.TypedArray.JSONEncodedArray: JSONEncodedArray.JSONEncodedArray,
    }
    (u_type, u) = tarr
    if u_type is TypedArray.TypedArray.NONE:
        return None

    TarType = type_map.get(u_type, None)
    if TarType is None:
        raise TypeError(f"FBS contains unknown data type: {u_type}")

    arr = TarType()
    arr.Init(u.Bytes, u.Pos)
    narr = arr.DataAsNumpy()
    if u_type == TypedArray.TypedArray.JSONEncodedArray:
        narr = json.loads(narr.tobytes().decode("utf-8"))
    return narr


def decode_matrix_fbs(fbs):
    """
    Given an FBS-encoded Matrix, return a Pandas DataFrame the contains the data and indices.
    """

    matrix = Matrix.Matrix.GetRootAsMatrix(fbs, 0)
    n_rows = matrix.NRows()
    n_cols = matrix.NCols()
    if n_rows == 0 or n_cols == 0:
        return pd.DataFrame()

    if matrix.RowIndexType() is not TypedArray.TypedArray.NONE:
        raise ValueError("row indexing not supported for FBS Matrix")

    columns_length = matrix.ColumnsLength()

    columns_index = deserialize_typed_array((matrix.ColIndexType(), matrix.ColIndex()))
    if columns_index is None:
        columns_index = range(0, n_cols)

    # sanity checks
    if len(columns_index) != n_cols or columns_length != n_cols:
        raise ValueError("FBS column count does not match number of columns in underlying matrix")

    columns_data = {}
    columns_type = {}
    for col_idx in range(0, columns_length):
        col = matrix.Columns(col_idx)
        tarr = (col.UType(), col.U())
        data = deserialize_typed_array(tarr)
        columns_data[columns_index[col_idx]] = data
        if len(data) != n_rows:
            raise ValueError("FBS column length does not match number of rows")
        if col.UType() is TypedArray.TypedArray.JSONEncodedArray:
            columns_type[columns_index[col_idx]] = "category"

    df = pd.DataFrame.from_dict(data=columns_data).astype(columns_type, copy=False)

    # more sanity checks
    if not df.columns.is_unique or len(df.columns) != n_cols:
        raise KeyError("FBS column indices are not unique")

    return df
