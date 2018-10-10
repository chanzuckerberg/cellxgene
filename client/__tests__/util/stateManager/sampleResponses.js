/* eslint no-bitwise: "off" */
import _ from "lodash";

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

const anAnnotationsObsResponse = {
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

const anAnnotationsVarResponse = {
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

const aLayoutResponse = {
  layout: {
    ndims: 2,
    coordinates: _()
      .range(nObs)
      .map(idx => [idx, Math.random(), Math.random()])
      .value()
  }
};

const aDataObsResponse = {
  var: [2, 4, 29],
  obs: _()
    .range(nObs)
    .map(idx => [idx, Math.random(), Math.random(), Math.random()])
    .value()
};

export {
  aLayoutResponse as layoutObs,
  aDataObsResponse as dataObs,
  anAnnotationsVarResponse as annotationsVar,
  anAnnotationsObsResponse as annotationsObs,
  aSchemaResponse as schema,
  aConfigResponse as config
};
