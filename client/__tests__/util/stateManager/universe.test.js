import * as Universe from "../../../src/util/stateManager/universe";
import * as Dataframe from "../../../src/util/dataframe";
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
        obsAnnotations: expect.any(Dataframe.Dataframe),
        varAnnotations: expect.any(Dataframe.Dataframe),
        obsLayout: expect.any(Dataframe.Dataframe),
        summary: expect.any(Object),
        varDataCache: expect.any(Object)
      })
    );

    expect(universe.obsAnnotations.dims).toEqual([
      nObs,
      REST.schema.schema.annotations.obs.length
    ]);
    expect(universe.obsLayout.dims).toEqual([nObs, 2]);
    expect(universe.obsLayout.colIndex.keys()).toEqual(["X", "Y"]);
    expect(universe.varAnnotations.dims).toEqual([
      nVar,
      REST.schema.schema.annotations.var.length
    ]);
  });
});
