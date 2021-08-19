import every from "lodash.every";
import map from "lodash.map";
import isNumber from "lodash.isnumber";
import zip from "lodash.zip";
import _ from "lodash";
import { flatbuffers } from "flatbuffers";
import { NetEncoding } from "../../../src/util/stateManager/matrix_generated";

/*
test data mocking REST 0.2 API responses.  Used in several tests.
*/
const nObs = 10;
const nVar = 32;
const field4Categories = [83, true, "foo", 2.222222];
const fieldDCategories = [99, false, "mumble", 3.1415];

const aConfigResponse = {
  config: {
    features: [
      { method: "POST", path: "/cluster/", available: false },
      { method: "POST", path: "/layout/", available: false },
      { method: "POST", path: "/diffexp/", available: false },
      { method: "POST", path: "/saveLocal/", available: false },
    ],
    displayNames: {
      engine: "the little engine that could",
      dataset: "all your zeros are mine",
    },
  },
};

const aSchemaResponse = {
  schema: {
    dataframe: {
      nObs,
      nVar,
      type: "float32",
    },
    annotations: {
      obs: {
        index: "name",
        columns: [
          { name: "name", type: "string" },
          { name: "field1", type: "int32" },
          { name: "field2", type: "float32" },
          { name: "field3", type: "boolean" },
          {
            name: "field4",
            type: "categorical",
            categories: field4Categories,
          },
        ],
      },
      var: {
        index: "name",
        columns: [
          { name: "name", type: "string" },
          { name: "fieldA", type: "int32" },
          { name: "fieldB", type: "float32" },
          { name: "fieldC", type: "boolean" },
          {
            name: "fieldD",
            type: "categorical",
            categories: fieldDCategories,
          },
        ],
      },
    },
    layout: {
      obs: [{ name: "umap", type: "float32", dims: ["umap_0", "umap_1"] }],
      var: [],
    },
  },
};

const anAnnotationsObsJSONResponse = {
  names: ["name", "field1", "field2", "field3", "field4"],
  // @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 0.
  data: _()
    .range(nObs)
    .map((idx) => [
      idx,
      `obs${idx}`,
      2 * idx,
      idx + 0.0133,
      // eslint-disable-next-line no-bitwise -- idx & 1 to check for odd numbers
      !!(idx & 1),
      field4Categories[idx % field4Categories.length],
    ])
    .value(),
};

const anAnnotationsVarJSONResponse = {
  names: ["fieldA", "fieldB", "fieldC", "fieldD", "name"],
  // @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 0.
  data: _()
    .range(nVar)
    .map((idx) => [
      idx,
      10 * idx,
      idx + 2.90143,
      // eslint-disable-next-line no-bitwise -- idx & 1 to check for odd numbers
      !!(idx & 1),
      fieldDCategories[idx % fieldDCategories.length],
      `var${idx}`,
    ])
    .value(),
};

function encodeTypedArray(builder: any, uType: any, uData: any) {
  const uTypeName = NetEncoding.TypedArray[uType];
  const ArrayType = NetEncoding[uTypeName];
  const dv = ArrayType.createDataVector(builder, uData);
  builder.startObject(1);
  builder.addFieldOffset(0, dv, 0);
  return builder.endObject();
}

function encodeMatrix(columns: any, colIndex = undefined) {
  /*
  IMPORTANT: this is not a general purpose encoder.  in particular,
  it doesn't correctly handle all column index types, nor does it
  handle all column typedarray types.

  encodeMatrixFBS in matrix.py is more general.  This is used only
  as a testing santity check (alt implementation).
  */
  // @ts-expect-error ts-migrate(2554) FIXME: Expected 0 arguments, but got 1.
  const utf8Encoder = new TextEncoder("utf-8");
  const builder = new flatbuffers.Builder(1024);
  const cols = map(columns, (carr) => {
    let uType;
    let tarr;
    if (every(carr, isNumber)) {
      uType = NetEncoding.TypedArray.Float32Array;
      tarr = encodeTypedArray(builder, uType, new Float32Array(carr));
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

  const encColumns = NetEncoding.Matrix.createColumnsVector(builder, cols);

  let encColIndex;
  if (colIndex) {
    encColIndex = encodeTypedArray(
      builder,
      NetEncoding.TypedArray.JSONEncodedArray,
      utf8Encoder.encode(JSON.stringify(colIndex))
    );
  }

  NetEncoding.Matrix.startMatrix(builder);
  NetEncoding.Matrix.addNRows(builder, columns[0].length);
  NetEncoding.Matrix.addNCols(builder, columns.length);
  NetEncoding.Matrix.addColumns(builder, encColumns);
  if (colIndex) {
    NetEncoding.Matrix.addColIndexType(
      builder,
      NetEncoding.TypedArray.JSONEncodedArray
    );
    NetEncoding.Matrix.addColIndex(builder, encColIndex);
  }
  const root = NetEncoding.Matrix.endMatrix(builder);
  builder.finish(root);
  return builder.asUint8Array();
}

const anAnnotationsObsFBSResponse = (() => {
  const columns = zip(...anAnnotationsObsJSONResponse.data).slice(1);
  // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'string[]' is not assignable to p... Remove this comment to see the full error message
  return encodeMatrix(columns, anAnnotationsObsJSONResponse.names);
})();

const anAnnotationsVarFBSResponse = (() => {
  const columns = zip(...anAnnotationsVarJSONResponse.data).slice(1);
  // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'string[]' is not assignable to p... Remove this comment to see the full error message
  return encodeMatrix(columns, anAnnotationsVarJSONResponse.names);
})();

const aLayoutFBSResponse = (() => {
  const coords = [
    new Float32Array(nObs).fill(Math.random()),
    new Float32Array(nObs).fill(Math.random()),
  ];
  // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'string[]' is not assignable to p... Remove this comment to see the full error message
  return encodeMatrix(coords, ["umap_0", "umap_1"]);
})();

const aDataObsResponse = {
  var: [2, 4, 29],
  // @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 0.
  obs: _()
    .range(nObs)
    .map((idx) => [idx, Math.random(), Math.random(), Math.random()])
    .value(),
};

export {
  aLayoutFBSResponse as layoutObs,
  aDataObsResponse as dataObs,
  anAnnotationsVarFBSResponse as annotationsVar,
  anAnnotationsObsFBSResponse as annotationsObs,
  aSchemaResponse as schema,
  aConfigResponse as config,
};
