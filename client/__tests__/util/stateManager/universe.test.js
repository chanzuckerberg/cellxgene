import _ from "lodash";
import * as Universe from "../../../src/util/stateManager/universe";
import * as REST from "./sampleResponses";

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

    expect(universe.obsAnnotations).toHaveLength(nObs);
    expect(_.keys(universe.obsNameToIndexMap)).toHaveLength(nObs);
    expect(universe.obsLayout.X).toHaveLength(nObs);
    expect(universe.obsLayout.Y).toHaveLength(nObs);
    expect(universe.varAnnotations).toHaveLength(nVar);
    expect(_.keys(universe.varNameToIndexMap)).toHaveLength(nVar);
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
