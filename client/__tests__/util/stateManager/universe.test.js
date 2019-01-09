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
