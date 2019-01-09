import flatbuffers
import numpy as np
from scipy import sparse
import pandas as pd

import server.app.util.fbs.NetEncoding.Column as Column
import server.app.util.fbs.NetEncoding.TypedArray as TypedArray
import server.app.util.fbs.NetEncoding.Matrix as Matrix


# Placeholder until recent enhancements to flatbuffers Python
# runtime are released, at which point we can use the default
# version.  This code is a port of the head.  See:
#
#       https://github.com/google/flatbuffers/pull/4829
#
def CreateNumpyVector(builder, x):
    """CreateNumpyVector writes a numpy array into the buffer."""

    if not isinstance(x, np.ndarray):
        raise TypeError("non-numpy-ndarray passed to CreateNumpyVector")

    if x.dtype.kind not in ['b', 'i', 'u', 'f']:
        raise TypeError("numpy-ndarray holds elements of unsupported datatype")

    if x.ndim > 1:
        raise TypeError("multidimensional-ndarray passed to CreateNumpyVector")

    builder.StartVector(x.itemsize, x.size, x.dtype.alignment)

    # Ensure little endian byte ordering
    if x.dtype.str[0] == "<":
        x_little_endian = x
    else:
        x_little_endian = x.byteswap(inplace=False)

    # Calculate total length
    len = int(x_little_endian.itemsize * x_little_endian.size)
    builder.head = int(builder.Head() - len)

    # tobytes ensures c_contiguous ordering
    builder.Bytes[builder.Head():builder.Head() + len] = x_little_endian.tobytes(order='C')

    return builder.EndVector(x.size)


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
    Serialize any of the various typed arrays, eg, Float32Array.   Specific
    means of serialization and type conversion are provided by type_info.
    """
    arr = source_array
    (array_type, as_type) = encoding_info(source_array)

    if isinstance(arr, pd.Index):
        arr = arr.to_series()

    # convert to a simple ndarray
    if as_type == 'json':
        as_json = arr.to_json(orient='records')
        arr = np.array(bytearray(as_json, 'utf-8'))
    else:
        if sparse.issparse(arr):
            arr = arr.toarray()
        elif isinstance(arr, pd.Series):
            arr = arr.get_values()
        if arr.dtype != as_type:
            arr = arr.astype(as_type)

    # serialize the ndarray into a vector
    if arr.ndim == 2 and arr.shape[0] == 1:
        arr = arr[0]
    vec = CreateNumpyVector(builder, arr)

    # serialize the typed array table
    builder.StartObject(1)
    builder.PrependUOffsetTRelativeSlot(0, vec, 0)
    array_value = builder.EndObject()
    return (array_type, array_value)


def column_encoding(arr):
    type_map = {
        # dtype:  ( array_type, as_type )
        np.float64: (TypedArray.TypedArray.Float32Array, np.float32),
        np.float32: (TypedArray.TypedArray.Float32Array, np.float32),
        np.float16: (TypedArray.TypedArray.Float32Array, np.float32),

        np.int8: (TypedArray.TypedArray.Int32Array, np.int32),
        np.int16: (TypedArray.TypedArray.Int32Array, np.int32),
        np.int32: (TypedArray.TypedArray.Int32Array, np.int32),
        np.int64: (TypedArray.TypedArray.Int32Array, np.int32),

        np.uint8: (TypedArray.TypedArray.Uint32Array, np.uint32),
        np.uint16: (TypedArray.TypedArray.Uint32Array, np.uint32),
        np.uint32: (TypedArray.TypedArray.Uint32Array, np.uint32),
        np.uint64: (TypedArray.TypedArray.Uint32Array, np.uint32)
    }
    type_map_default = (TypedArray.TypedArray.JSONEncodedArray, 'json')
    return type_map.get(arr.dtype.type, type_map_default)


def index_encoding(arr):
    type_map = {
        # dtype:  ( array_type, as_type )
        np.int32: (TypedArray.TypedArray.Int32Array, np.int32),
        np.int64: (TypedArray.TypedArray.Int32Array, np.int32),

        np.uint32: (TypedArray.TypedArray.Uint32Array, np.uint32),
        np.uint64: (TypedArray.TypedArray.Uint32Array, np.uint32)
    }
    type_map_default = (TypedArray.TypedArray.JSONEncodedArray, 'json')
    return type_map.get(arr.dtype.type, type_map_default)


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
    guess = (guess + 0x400) & (~0x3ff)
    return guess


def encode_matrix_fbs(matrix, row_idx=None, col_idx=None):
    """
    Given a 2D DataFrame, ndarray or sparse equivalent, create and return a
    Matrix flatbuffer.

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
    builder = flatbuffers.Builder(guess_at_mem_needed(matrix))

    if isinstance(matrix, pd.DataFrame):
        matrix_columns = reversed(tuple(matrix[name] for name in matrix))
    else:
        matrix_columns = reversed(tuple(c for c in matrix.T))

    columns = []
    # for idx in reversed(np.arange(n_cols)):
    for c in matrix_columns:
        # serialize the typed array
        typed_arr = serialize_typed_array(builder, c, column_encoding)

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
