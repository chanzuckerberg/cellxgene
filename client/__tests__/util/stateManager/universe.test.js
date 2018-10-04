import _ from "lodash";
import * as Universe from "../../../src/util/stateManager/universe";

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

  const aSchemaResponse = {
    schema: {

    };
  };

  const aConfigResponse = {
    config: {
      
    }
  }

  test("", () => {
    /*
    */
  });
});

describe("convertExpressionRESTv02ToObject", () => {
  /*
  test convertExpressionRESTv02ToObject
  */
});
