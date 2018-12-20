/* eslint no-bitwise: "off" */
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
      { method: "POST", path: "/saveLocal/", available: false }
    ],
    displayNames: {
      engine: "the little engine that could",
      dataset: "all your zeros are mine"
    }
  }
};

const aSchemaResponse = {
  schema: {
    dataframe: {
      nObs,
      nVar,
      type: "float32"
    },
    annotations: {
      obs: [
        { name: "name", type: "string" },
        { name: "field1", type: "int32" },
        { name: "field2", type: "float32" },
        { name: "field3", type: "boolean" },
        {
          name: "field4",
          type: "categorical",
          categories: field4Categories
        }
      ],
      var: [
        { name: "name", type: "string" },
        { name: "fieldA", type: "int32" },
        { name: "fieldB", type: "float32" },
        { name: "fieldC", type: "boolean" },
        {
          name: "fieldD",
          type: "categorical",
          categories: fieldDCategories
        }
      ]
    }
  }
};

const anAnnotationsObsJSONResponse = {
  names: ["name", "field1", "field2", "field3", "field4"],
  data: _()
    .range(nObs)
    .map(idx => [
      idx,
      `obs${idx}`,
      2 * idx,
      idx + 0.0133,
      !!(idx & 1),
      field4Categories[idx % field4Categories.length]
    ])
    .value()
};

const anAnnotationsVarJSONResponse = {
  names: ["fieldA", "fieldB", "fieldC", "fieldD", "name"],
  data: _()
    .range(nVar)
    .map(idx => [
      idx,
      10 * idx,
      idx + 2.90143,
      !!(idx & 1),
      fieldDCategories[idx % fieldDCategories.length],
      `var${idx}`
    ])
    .value()
};

function encodeTypedArray(builder, uType, uData) {
  const uTypeName = NetEncoding.TypedArray[uType];
  const ArrayType = NetEncoding[uTypeName];
  const dv = ArrayType.createDataVector(builder, uData);
  builder.startObject(1);
  builder.addFieldOffset(0, dv, 0);
  return builder.endObject();
}

function encodeDataFrame(columns, colIndex = undefined) {
  const utf8Encoder = new TextEncoder("utf-8");
  const builder = new flatbuffers.Builder(1024);
  const cols = _.map(columns, carr => {
    let uType;
    let tarr;
    if (_.every(carr, _.isNumber)) {
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

  const encColumns = NetEncoding.DataFrame.createColumnsVector(builder, cols);

  let encColIndex;
  if (colIndex) {
    encColIndex = encodeTypedArray(
      builder,
      NetEncoding.TypedArray.JSONEncodedArray,
      utf8Encoder.encode(JSON.stringify(colIndex))
    );
  }

  NetEncoding.DataFrame.startDataFrame(builder);
  NetEncoding.DataFrame.addNRows(builder, columns[0].length);
  NetEncoding.DataFrame.addNCols(builder, columns.length);
  NetEncoding.DataFrame.addColumns(builder, encColumns);
  if (colIndex) {
    NetEncoding.DataFrame.addColIndexType(
      builder,
      NetEncoding.TypedArray.JSONEncodedArray
    );
    NetEncoding.DataFrame.addColIndex(builder, encColIndex);
  }
  const root = NetEncoding.DataFrame.endDataFrame(builder);
  builder.finish(root);
  return builder.asUint8Array();
}

const anAnnotationsObsFBSResponse = (() => {
  const columns = _.zip(...anAnnotationsObsJSONResponse.data).slice(1);
  return encodeDataFrame(columns, anAnnotationsObsJSONResponse.names);
})();

const anAnnotationsVarFBSResponse = (() => {
  const columns = _.zip(...anAnnotationsVarJSONResponse.data).slice(1);
  return encodeDataFrame(columns, anAnnotationsVarJSONResponse.names);
})();

const aLayoutJSONResponse = {
  layout: {
    ndims: 2,
    coordinates: _()
      .range(nObs)
      .map(idx => [idx, Math.random(), Math.random()])
      .value()
  }
};

const aLayoutFBSResponse = (() => {
  const coords = [
    new Float32Array(nObs).fill(Math.random()),
    new Float32Array(nObs).fill(Math.random())
  ];
  const builder = new flatbuffers.Builder(1024);

  const cols = _.map(coords, carr => {
    const cdv = NetEncoding.Float32Array.createDataVector(builder, carr);
    NetEncoding.Float32Array.startFloat32Array(builder);
    NetEncoding.Float32Array.addData(builder, cdv);
    const floatArr = NetEncoding.Float32Array.endFloat32Array(builder);

    NetEncoding.Column.startColumn(builder);
    NetEncoding.Column.addUType(builder, NetEncoding.TypedArray.Float32Array);
    NetEncoding.Column.addU(builder, floatArr);
    return NetEncoding.Column.endColumn(builder);
  });

  const columns = NetEncoding.DataFrame.createColumnsVector(builder, cols);

  NetEncoding.DataFrame.startDataFrame(builder);
  NetEncoding.DataFrame.addNRows(builder, nObs);
  NetEncoding.DataFrame.addNCols(builder, nVar);
  NetEncoding.DataFrame.addColumns(builder, columns);
  const matrix = NetEncoding.DataFrame.endDataFrame(builder);
  builder.finish(matrix);
  return builder.asUint8Array();
})();

const aDataObsResponse = {
  var: [2, 4, 29],
  obs: _()
    .range(nObs)
    .map(idx => [idx, Math.random(), Math.random(), Math.random()])
    .value()
};

export {
  aLayoutFBSResponse as layoutObs,
  aDataObsResponse as dataObs,
  anAnnotationsVarFBSResponse as annotationsVar,
  anAnnotationsObsFBSResponse as annotationsObs,
  aSchemaResponse as schema,
  aConfigResponse as config
};
