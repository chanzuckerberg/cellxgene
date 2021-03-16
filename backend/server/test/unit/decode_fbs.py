"""
Code to decode, for testing purposes, the flatbuffer encoded blobs.

This code will need to be updated if fbs/matrix.fbs changes. For more information, see fbs/matrix.fbs and
server/data_common/fbs/
"""

import backend.common_utils.fbs.NetEncoding.Matrix as Matrix
from backend.common_utils.fbs.matrix import deserialize_typed_array


def decode_matrix_FBS(buf):
    """
    Given a FBS Matrix, return an decoded Python dict containing same info in native format.
    NOTE / TODO: row_idx not currently implemented
    """
    df = Matrix.Matrix.GetRootAsMatrix(buf, 0)
    n_rows = df.NRows()
    n_cols = df.NCols()

    columns_length = df.ColumnsLength()

    decoded_columns = []
    for col_idx in range(0, columns_length):
        col = df.Columns(col_idx)
        tarr = (col.UType(), col.U())
        decoded_columns.append(deserialize_typed_array(tarr))

    cidx = deserialize_typed_array((df.ColIndexType(), df.ColIndex()))

    return {"n_rows": n_rows, "n_cols": n_cols, "columns": decoded_columns, "col_idx": cidx, "row_idx": None}
