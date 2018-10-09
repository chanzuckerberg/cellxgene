import _ from "lodash";
import * as Universe from "../../../src/util/stateManager/universe";
import * as REST from "./sampleResponses";

// /*
// test data mocking REST 0.2 API responses
// */
// const nObs = 10;
// const nVar = 32;
// const field4Categories = [83, true, "foo", 2.222222];
// const fieldDCategories = [99, false, "mumble", 3.1415];
//
// const aConfigResponse = {
//   config: {
//     features: [
//       { method: "POST", path: "/cluster/", available: false },
//       { method: "POST", path: "/layout/", available: false },
//       { method: "POST", path: "/diffexp/", available: false },
//       { method: "POST", path: "/saveLocal/", available: false }
//     ],
//     displayNames: {
//       engine: "the little engine that could",
//       dataset: "all your zeros are mine"
//     }
//   }
// };
//
// const aSchemaResponse = {
//   schema: {
//     dataframe: {
//       nObs,
//       nVar,
//       type: "float32"
//     },
//     annotations: {
//       obs: [
//         { name: "name", type: "string" },
//         { name: "field1", type: "int32" },
//         { name: "field2", type: "float32" },
//         { name: "field3", type: "boolean" },
//         {
//           name: "field4",
//           type: "categorical",
//           categories: field4Categories
//         }
//       ],
//       var: [
//         { name: "name", type: "string" },
//         { name: "fieldA", type: "int32" },
//         { name: "fieldB", type: "float32" },
//         { name: "fieldC", type: "boolean" },
//         {
//           name: "fieldD",
//           type: "categorical",
//           categories: fieldDCategories
//         }
//       ]
//     }
//   }
// };
//
// const anAnnotationsObsResponse = {
//   names: ["name", "field1", "field2", "field3", "field4"],
//   data: _()
//     .range(nObs)
//     .map(idx => [
//       idx,
//       `obs${idx}`,
//       2 * idx,
//       idx + 0.0133,
//       idx & 1,
//       field4Categories[idx % field4Categories.length]
//     ])
//     .value()
// };
//
// const anAnnotationsVarResponse = {
//   names: ["fieldA", "fieldB", "fieldC", "fieldD", "name"],
//   data: _()
//     .range(nVar)
//     .map(idx => [
//       idx,
//       10 * idx,
//       idx + 2.90143,
//       idx & 1,
//       fieldDCategories[idx % fieldDCategories.length],
//       `var${idx}`
//     ])
//     .value()
// };
//
// const aLayoutResponse = {
//   layout: {
//     ndims: 2,
//     coordinates: _()
//       .range(nObs)
//       .map(idx => [idx, Math.random(), Math.random()])
//       .value()
//   }
// };
//
// const aDataObsResponse = {
//   var: [2, 4, 29],
//   obs: _()
//     .range(nObs)
//     .map(idx => [idx, Math.random(), Math.random(), Math.random()])
//     .value()
// };

describe("createUniverseFromRestV02Response", () => {
  /*
  test createUniverseFromRestV02Response - this function converts
  a set of REST 0.2 responses into a "new" Universe.

  createUniverseFromRestV02Response(
    configResponse,
    schemaResponse,
    annotationsObsResponse,
    annotationsVarResponse,
    layoutObsResponse
  ) --> Universe

  where:
  configResponse: GET /.../config
  schemaResponse: GET /.../schema
  annotationsObsResponse: GET /.../annotations/obs
  annotationsVarResponse: GET /.../annotations/var
  layoutObsResponse: GET /.../layout/obs

  See spec in docs/REST_API.md.
  */

  test("create from test data", () => {
    /*
    create a universe from sample data nad validate its shape & contents
    */
    const { nObs, nVar } = REST.schema.schema.dataframe;

    const universe = Universe.createUniverseFromRestV02Response(
      REST.config,
      REST.schema,
      REST.annotationsObs,
      REST.annotationsVar,
      REST.layoutObs
    );

    expect(universe).toBeDefined();
    expect(universe).toMatchObject(
      expect.objectContaining({
        api: "0.2",
        nObs,
        nVar,
        schema: REST.schema.schema,
        obsAnnotations: expect.any(Array),
        varAnnotations: expect.any(Array),
        obsNameToIndexMap: expect.any(Object),
        varNameToIndexMap: expect.any(Object),
        obsLayout: expect.objectContaining({
          X: expect.any(Float32Array),
          Y: expect.any(Float32Array)
        }),
        varDataCache: expect.any(Object)
      })
    );

    expect(_.size(universe.obsAnnotations)).toBe(nObs);
    expect(_.size(universe.obsNameToIndexMap)).toBe(nObs);
    expect(_.size(universe.obsLayout.X)).toBe(nObs);
    expect(_.size(universe.obsLayout.Y)).toBe(nObs);
    expect(_.size(universe.varAnnotations)).toBe(nVar);
    expect(_.size(universe.varNameToIndexMap)).toBe(nVar);
  });
});

describe("convertExpressionRESTv02ToObject", () => {
  /*
  test convertExpressionRESTv02ToObject

  convertExpressionRESTv02ToObject(
    universe,
    response) --> { geneName: Float32Array, geneName: Float32Array, ... }

  reponse is a /data/obs response:
  {
    var: [ varIndices fetched ],
    obs: [
      [ obsIndex, evalue, ... ],
      ...
    ]
  }
  */
  test("create from response data", () => {
    const universe = Universe.createUniverseFromRestV02Response(
      REST.config,
      REST.schema,
      REST.annotationsObs,
      REST.annotationsVar,
      REST.layoutObs
    );
    const expression = Universe.convertExpressionRESTv02ToObject(
      universe,
      REST.dataObs
    );

    /* Check that the expected keys are present */
    const expectedGeneNames = _.map(
      REST.dataObs.var,
      v => REST.annotationsVar.data[v][5]
    );
    expect(Object.keys(expression)).toEqual(
      expect.arrayContaining(expectedGeneNames)
    );

    const expectedExpressionValues = _.map(
      _.unzip(REST.dataObs.obs),
      a => new Float32Array(a)
    );

    _.forEach(REST.dataObs.var, (varIdx, idx) => {
      const varName = universe.varAnnotations[varIdx].name;
      expect(varName).toBeDefined();
      expect(varIdx).toBe(universe.varNameToIndexMap[varName]);
      expect(expression[varName]).toEqual(expectedExpressionValues[idx + 1]);
    });
  });
});
