import { flatbuffers } from "flatbuffers";
import { NetEncoding } from "./matrix_generated";
import { isTypedArray, isFpTypedArray } from "../typeHelpers";
import {
  Dataframe,
  IdentityInt32Index,
  DenseInt32Index,
  KeyIndex,
} from "../dataframe";

const utf8Decoder = new TextDecoder("utf-8");

/*
Matrix flatbuffer decoding support.   See fbs/matrix.fbs
*/

/*
Decode NetEncoding.TypedArray
*/
function decodeTypedArray(uType, uValF, inplace = false) {
  if (uType === NetEncoding.TypedArray.NONE) {
    return null;
  }

  // Convert to a JS class that supports this type
  const TypeClass = NetEncoding[NetEncoding.TypedArray[uType]];
  // Create a TypedArray that references the underlying buffer
  let arr = uValF(new TypeClass()).dataArray();
  if (uType === NetEncoding.TypedArray.JSONEncodedArray) {
    const json = utf8Decoder.decode(arr);
    arr = JSON.parse(json);
  } else if (!inplace) {
    /* force copy to release underlying FBS buffer */
    arr = new arr.constructor(arr);
  }
  return arr;
}

/*
Parameter: Uint8Array or ArrayBuffer containing raw flatbuffer Matrix
Returns: object containing decoded Matrix:
{
  nRows: num,
  nCols: num,
  columns: [
    each column, which will be a TypedArray or Array
  ]
  colIdx: []|null
}
*/
export function decodeMatrixFBS(arrayBuffer, inplace = false) {
  const bb = new flatbuffers.ByteBuffer(new Uint8Array(arrayBuffer));
  const matrix = NetEncoding.Matrix.getRootAsMatrix(bb);

  const nRows = matrix.nRows();
  const nCols = matrix.nCols();

  /* decode columns */
  const columnsLength = matrix.columnsLength();
  const columns = Array(columnsLength).fill(null);
  for (let c = 0; c < columnsLength; c += 1) {
    const col = matrix.columns(c);
    columns[c] = decodeTypedArray(col.uType(), col.u.bind(col), inplace);
  }

  /* decode col_idx */
  const colIdx = decodeTypedArray(
    matrix.colIndexType(),
    matrix.colIndex.bind(matrix),
    inplace
  );

  return {
    nRows,
    nCols,
    columns,
    colIdx,
    rowIdx: null,
  };
}

function encodeTypedArray(builder, uType, uData) {
  const uTypeName = NetEncoding.TypedArray[uType];
  const ArrayType = NetEncoding[uTypeName];
  const dv = ArrayType.createDataVector(builder, uData);
  builder.startObject(1);
  builder.addFieldOffset(0, dv, 0);
  return builder.endObject();
}

export function encodeMatrixFBS(df) {
  /*
  encode the dataframe as an FBS Matrix
  */

  /* row indexing not supported currently */
  if (df.rowIndex.constructor !== IdentityInt32Index) {
    throw new Error("FBS does not support row index encoding at this time");
  }

  const shape = df.dims;
  const utf8Encoder = new TextEncoder("utf-8");
  const builder = new flatbuffers.Builder(1024);

  let encColIndex;
  let encColIndexUType;
  let encColumns;

  if (shape[0] > 0 && shape[1] > 0) {
    const columns = df.columns().map((col) => col.asArray());

    const cols = columns.map((carr) => {
      let uType;
      let tarr;
      if (isTypedArray(carr)) {
        uType = NetEncoding.TypedArray[carr.constructor.name];
        tarr = encodeTypedArray(builder, uType, carr);
      } else {
        uType = NetEncoding.TypedArray.JSONEncodedArray;
        const json = JSON.stringify(carr);
        const jsonUTF8 = utf8Encoder.encode(json);
        tarr = encodeTypedArray(builder, uType, jsonUTF8);
      }
      NetEncoding.Column.startColumn(builder);
      NetEncoding.Column.addUType(builder, uType);
      NetEncoding.Column.addU(builder, tarr);
      return NetEncoding.Column.endColumn(builder);
    });

    encColumns = NetEncoding.Matrix.createColumnsVector(builder, cols);

    if (df.colIndex && shape[1] > 0) {
      const colIndexType = df.colIndex.constructor;
      if (colIndexType === IdentityInt32Index) {
        encColIndex = undefined;
      } else if (colIndexType === DenseInt32Index) {
        encColIndexUType = NetEncoding.TypedArray.Int32Array;
        encColIndex = encodeTypedArray(
          builder,
          encColIndexUType,
          df.colIndex.labels()
        );
      } else if (colIndexType === KeyIndex) {
        encColIndexUType = NetEncoding.TypedArray.JSONEncodedArray;
        encColIndex = encodeTypedArray(
          builder,
          encColIndexUType,
          utf8Encoder.encode(JSON.stringify(df.colIndex.labels()))
        );
      } else {
        throw new Error("Index type FBS encoding unsupported");
      }
    }
  }

  NetEncoding.Matrix.startMatrix(builder);
  NetEncoding.Matrix.addNRows(builder, shape[0]);
  NetEncoding.Matrix.addNCols(builder, shape[1]);
  if (encColumns) {
    NetEncoding.Matrix.addColumns(builder, encColumns);
  }
  if (encColIndexUType) {
    NetEncoding.Matrix.addColIndexType(builder, encColIndexUType);
    NetEncoding.Matrix.addColIndex(builder, encColIndex);
  }
  const root = NetEncoding.Matrix.endMatrix(builder);
  builder.finish(root);
  return builder.asUint8Array();
}

function promoteTypedArray(o) {
  /*
  Decide what internal data type to use for the data returned from
  the server.

  TODO - future optimization: not all int32/uint32 data series require
  promotion to float64.  We COULD simply look at the data to decide.
  */
  if (isFpTypedArray(o) || Array.isArray(o)) return o;

  let TyepdArrayCtor;
  switch (o.constructor) {
    case Int8Array:
    case Uint8Array:
    case Uint8ClampedArray:
    case Int16Array:
    case Uint16Array:
      TyepdArrayCtor = Float32Array;
      break;

    case Int32Array:
    case Uint32Array:
      TyepdArrayCtor = Float64Array;
      break;

    default:
      throw new Error("Unexpected data type returned from server.");
  }
  if (o.constructor === TyepdArrayCtor) return o;
  return new TyepdArrayCtor(o);
}

export function matrixFBSToDataframe(arrayBuffers) {
  /*
  Convert array of Matrix FBS to a Dataframe.

  The application has strong assumptions that all scalar data will be
  stored as a float32 or float64 (regardless of underlying data types).
  For example, clipping of value ranges (eg, user-selected percentiles)
  depends on the ability to use NaN in any numeric type.

  All float data from the server is left as is.  All non-float is promoted
  to an appropriate float.
  */
  if (!Array.isArray(arrayBuffers)) {
    arrayBuffers = [arrayBuffers];
  }
  if (arrayBuffers.length === 0) {
    return Dataframe.Dataframe.empty();
  }

  const fbs = arrayBuffers.map((ab) => decodeMatrixFBS(ab, true)); // leave in place
  /* check that all FBS have same row dimensionality */
  const { nRows } = fbs[0];
  fbs.forEach((b) => {
    if (b.nRows !== nRows)
      throw new Error("FBS with inconsistent dimensionality");
  });
  const columns = fbs
    .map((fb) =>
      fb.columns.map((c) => {
        if (isFpTypedArray(c) || Array.isArray(c)) return c;
        return promoteTypedArray(c);
      })
    )
    .flat();
  // colIdx may be TypedArray or Array
  const colIdx = fbs
    .map((b) => (Array.isArray(b.colIdx) ? b.colIdx : Array.from(b.colIdx)))
    .flat();
  const nCols = columns.length;

  const df = new Dataframe([nRows, nCols], columns, null, new KeyIndex(colIdx));
  return df;
}
