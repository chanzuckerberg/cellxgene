/*
test FBS encode/decode API
*/
import { Dataframe, KeyIndex } from "../../../src/util/dataframe";
import {
  decodeMatrixFBS,
  encodeMatrixFBS,
} from "../../../src/util/stateManager/matrix";

describe("encode/decode", () => {
  test("round trip", () => {
    const columns = [
      ["red", "green", "blue"],
      new Int32Array(3).fill(0),
      new Uint32Array(3).fill(1),
      new Float32Array(3).fill(2),
    ];

    const dfNoColIdx = new Dataframe([3, 4], columns);
    const dfA = decodeMatrixFBS(encodeMatrixFBS(dfNoColIdx));
    expect([dfA.nRows, dfA.nCols]).toEqual(dfNoColIdx.dims);
    expect(dfA.colIdx).toBeNull();
    expect(dfA.rowIdx).toBeNull();
    expect(dfA.columns).toEqual(columns);

    const colIndex = new KeyIndex(["a", "b", "c", "d"]);
    const dfWithColIdx = new Dataframe([3, 4], columns, null, colIndex);
    const dfB = decodeMatrixFBS(encodeMatrixFBS(dfWithColIdx));
    expect([dfB.nRows, dfB.nCols]).toEqual(dfWithColIdx.dims);
    expect(dfB.colIdx).toEqual(colIndex.labels());
    expect(dfB.rowIdx).toBeNull();
    expect(dfB.columns).toEqual(columns);
  });
});
