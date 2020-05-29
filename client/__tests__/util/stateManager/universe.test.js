import * as Universe from "../../../src/util/stateManager/universe";
import { matrixFBSToDataframe } from "../../../src/util/stateManager/matrix";
import * as Dataframe from "../../../src/util/dataframe";
import * as REST from "./sampleResponses";

describe("createUniverseFromResponse", () => {
  /*
  test createUniverseFromResponse - this function converts
  a set of REST 0.2 responses into a "new" Universe.

  createUniverseFromResponse(
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
    let universe = Universe.createUniverseFromResponse(
      REST.config,
      REST.schema
    );
    expect(universe).toBeDefined();
    expect(universe).toMatchObject(
      expect.objectContaining({
        nObs,
        nVar,
        schema: REST.schema.schema,
        obsAnnotations: expect.any(Dataframe.Dataframe),
        varAnnotations: expect.any(Dataframe.Dataframe),
        obsLayout: expect.any(Dataframe.Dataframe),
        varData: expect.any(Dataframe.Dataframe),
      })
    );

    universe = {
      ...universe,
      ...Universe.addObsAnnotations(
        universe,
        matrixFBSToDataframe(REST.annotationsObs)
      ),
      ...Universe.addVarAnnotations(
        universe,
        matrixFBSToDataframe(REST.annotationsVar)
      ),
      ...Universe.addObsLayout(universe, matrixFBSToDataframe(REST.layoutObs)),
    };

    expect(universe).toMatchObject(
      expect.objectContaining({
        nObs,
        nVar,
        schema: REST.schema.schema,
        obsAnnotations: expect.any(Dataframe.Dataframe),
        varAnnotations: expect.any(Dataframe.Dataframe),
        obsLayout: expect.any(Dataframe.Dataframe),
        varData: expect.any(Dataframe.Dataframe),
      })
    );

    expect(universe.obsAnnotations.dims).toEqual([
      nObs,
      REST.schema.schema.annotations.obs.columns.length,
    ]);
    expect(universe.obsLayout.dims).toEqual([nObs, 2]);
    expect(universe.obsLayout.colIndex.labels()).toEqual(
      universe.schema.layout.obs[0].dims
    );
    expect(universe.varAnnotations.dims).toEqual([
      nVar,
      REST.schema.schema.annotations.var.columns.length,
    ]);
    expect(universe.varData.isEmpty()).toBeTruthy();
  });
});
