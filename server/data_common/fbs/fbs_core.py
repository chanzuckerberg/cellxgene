"""
Code to decode the flatbuffer encoded blobs in a python dict

This code may need to be updated if fbs/matrix.py changes.
"""

import json

import server.data_common.fbs.NetEncoding.Matrix as Matrix
import server.data_common.fbs.NetEncoding.Float32Array as Float32Array
import server.data_common.fbs.NetEncoding.Float64Array as Float64Array
import server.data_common.fbs.NetEncoding.Int32Array as Int32Array
import server.data_common.fbs.NetEncoding.JSONEncodedArray as JSONEncodedArray
import server.data_common.fbs.NetEncoding.TypedArray as TypedArray
import server.data_common.fbs.NetEncoding.Uint32Array as Uint32Array


def decode_matrix_fbs_to_dict(buf):
    """
    Given a FBS Matrix, return an decoded Python dict containing same info in native format.
    This function is mainly used for testing purposed.
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
        narr = json.loads(narr.tostring().decode("utf-8"))
    return narr
