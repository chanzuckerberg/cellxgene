"""
Code to decode, for testing purposes, the flatbuffer encoded blobs.
This code will need to be updated if fbs/matrix.fbs changes.

For more information, see fbs/matrix.fbs and server/data_common/fbs/
"""
import json

import server.data_common.fbs.NetEncoding.TypedArray as TypedArray
import server.data_common.fbs.NetEncoding.Matrix as Matrix
import server.data_common.fbs.NetEncoding.Int32Array as Int32Array
import server.data_common.fbs.NetEncoding.Uint32Array as Uint32Array
import server.data_common.fbs.NetEncoding.Float32Array as Float32Array
import server.data_common.fbs.NetEncoding.Float64Array as Float64Array
import server.data_common.fbs.NetEncoding.JSONEncodedArray as JSONEncodedArray


def decode_typed_array(tarr):
    type_map = {
        TypedArray.TypedArray.Uint32Array: Uint32Array.Uint32Array,
        TypedArray.TypedArray.Int32Array: Int32Array.Int32Array,
        TypedArray.TypedArray.Float32Array: Float32Array.Float32Array,
        TypedArray.TypedArray.Float64Array: Float64Array.Float64Array,
        TypedArray.TypedArray.JSONEncodedArray: JSONEncodedArray.JSONEncodedArray,
    }
    (u_type, u) = tarr
    if u_type == TypedArray.TypedArray.NONE:
        return None

    TarType = type_map.get(u_type, None)
    assert TarType is not None

    arr = TarType()
    arr.Init(u.Bytes, u.Pos)
    narr = arr.DataAsNumpy()
    if u_type == TypedArray.TypedArray.JSONEncodedArray:
        narr = json.loads(narr.tostring().decode("utf-8"))
    return narr


def decode_matrix_FBS(buf):
    """
    Given a FBS Matrix, return an decoded Python dict containing
    same info in native format.

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
        decoded_columns.append(decode_typed_array(tarr))

    cidx = decode_typed_array((df.ColIndexType(), df.ColIndex()))

    return {"n_rows": n_rows, "n_cols": n_cols, "columns": decoded_columns, "col_idx": cidx, "row_idx": None}
