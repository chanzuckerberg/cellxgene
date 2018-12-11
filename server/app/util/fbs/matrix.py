import flatbuffers
import numpy as np
from scipy import sparse

import server.app.util.fbs.NetEncoding.Column as Column
import server.app.util.fbs.NetEncoding.ColumnUnion as ColumnUnion
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
        x_lend = x
    else:
        x_lend = x.byteswap(inplace=False)

    # Calculate total length
    len = int(x_lend.itemsize * x_lend.size)
    builder.head = int(builder.Head() - len)

    # tobytes ensures c_contiguous ordering
    builder.Bytes[builder.Head():builder.Head()+len] = x_lend.tobytes(order='C')

    return builder.EndVector(x.size)


def create_matrix_flatbuffer(matrix, row_idx=None, col_idx=None):
    """
    Given an array of ndarrays, create and return a Matrix flatbuffer.

    :param matrix: 2D ndarray or sparse equivalent
    :param row_idx: numeric index for row dimension
    :param col_idx: numeric index for col dimension

    Currently indices are unsupported and must be None
    """

    if row_idx is not None:
        raise ValueError("row indexing not supproted for FBS Matrix")

    n_rows = matrix.shape[1]
    n_cols = matrix.shape[0]
    dtype = matrix.dtype

    # estimate size needed, so we don't unnecessarily realloc
    guess_at_mem = (n_rows * n_cols * dtype.itemsize) + 1024
    guess_at_mem = (guess_at_mem + 0x400) & (~0x3ff)  # round up to nearest 1024 bytes
    builder = flatbuffers.Builder(guess_at_mem)

    type_map = {
        # dtype:  ( array_type )
        np.float64: (ColumnUnion.ColumnUnion.Float64Array),
        np.float32: (ColumnUnion.ColumnUnion.Float32Array),
        np.int32: (ColumnUnion.ColumnUnion.Int32Array),
        np.uint32: (ColumnUnion.ColumnUnion.Uint32Array)
    }
    assert dtype.type in type_map
    array_type = type_map[dtype.type]

    columns = []
    for idx in reversed(np.arange(matrix.shape[0])):
        # serialize the underlying vector
        mat = matrix[idx]
        if sparse.issparse(mat):
            mat = mat.toarray()
        if mat.ndim == 2 and mat.shape[0] == 1:
            mat = mat[0]
        vec = CreateNumpyVector(builder, mat)

        # serialize the typed array
        builder.StartObject(1)
        builder.PrependUOffsetTRelativeSlot(0, vec, 0)
        array_value = builder.EndObject()

        # serialize the Column
        Column.ColumnStart(builder)
        Column.ColumnAddUType(builder, array_type)
        Column.ColumnAddU(builder, array_value)
        column = Column.ColumnEnd(builder)
        columns.append(column)

    # Serialize Matrix.columns[]
    Matrix.MatrixStartColumnsVector(builder, n_cols)
    for c in columns:
        builder.PrependUOffsetTRelative(c)
    mcv = builder.EndVector(n_cols)

    # serialize the colIndex if available
    cidx = None
    if col_idx is not None:
        cidx = CreateNumpyVector(builder, col_idx.astype('uint32'))

    # Serialize Matrix
    Matrix.MatrixStart(builder)
    Matrix.MatrixAddNRows(builder, n_rows)
    Matrix.MatrixAddNCols(builder, n_cols)
    Matrix.MatrixAddColumns(builder, mcv)
    if cidx is not None:
        Matrix.MatrixAddColIndex(builder, cidx)
    matrix = Matrix.MatrixEnd(builder)

    builder.Finish(matrix)
    return builder.Output()
