import { flatbuffers } from "flatbuffers";
import { NetEncoding } from "./matrix_generated";

const utf8Decoder = new TextDecoder("utf-8");

/*
DataFrame flatbuffer decoding support.   See fbs/matrix.fbs
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
Returns: object containing decoded DataFrame:
{
  nRows: num,
  nCols: num,
  columns: [
    each column, which will be a TypedArray or Array
  ]
  colIdx: []|null
}
*/
function decodeMatrixFBS(arrayBuffer, inplace = false) {
  const bb = new flatbuffers.ByteBuffer(new Uint8Array(arrayBuffer));
  const df = NetEncoding.DataFrame.getRootAsDataFrame(bb);

  const nRows = df.nRows();
  const nCols = df.nCols();

  /* decode columns */
  const columnsLength = df.columnsLength();
  const columns = Array(columnsLength).fill(null);
  for (let c = 0; c < columnsLength; c += 1) {
    const col = df.columns(c);
    columns[c] = decodeTypedArray(col.uType(), col.u.bind(col), inplace);
  }

  /* decode col_idx */
  const colIdx = decodeTypedArray(
    df.colIndexType(),
    df.colIndex.bind(df),
    inplace
  );

  return {
    nRows,
    nCols,
    columns,
    colIdx,
    rowIdx: null
  };
}

export default decodeMatrixFBS;
