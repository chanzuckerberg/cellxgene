import cloneDeep from "lodash.clonedeep";

import { NumberArray } from "../../src/common/types/arraytypes";
import calcCentroid from "../../src/util/centroid";
import quantile from "../../src/util/quantile";
import { matrixFBSToDataframe } from "../../src/util/stateManager/matrix";
import * as REST from "./stateManager/sampleResponses";
import { indexEntireSchema } from "../../src/util/stateManager/schemaHelpers";
import { normalizeWritableCategoricalSchema } from "../../src/annoMatrix/normalize";
import { Dataframe } from "../../src/util/dataframe";

describe("centroid", () => {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  let schema: any;
  let obsAnnotations: Dataframe;
  let obsLayout: Dataframe;

  beforeAll(() => {
    schema = indexEntireSchema(cloneDeep(REST.schema.schema));
    obsAnnotations = matrixFBSToDataframe(REST.annotationsObs);
    obsLayout = matrixFBSToDataframe(REST.layoutObs);

    normalizeWritableCategoricalSchema(
      schema.annotations.obsByName.field3,
      obsAnnotations.col("field3")
    );
  });

  test("field4 (categorical obsAnnotation)", () => {
    const centroidResult = calcCentroid(
      schema,
      "field4",
      obsAnnotations,
      { current: "umap", currentDimNames: ["umap_0", "umap_1"] },
      obsLayout
    );

    // Check to see that a centroid has been calculated for every categorical value
    const keysAsArray = Array.from(centroidResult.keys());
    expect(keysAsArray).toEqual(
      expect.arrayContaining([83, true, "foo", 2.222222])
    );

    // This expected result assumes that all cells belong in all categorical values inside of sample response
    const expectedResult = [
      quantile([0.5], obsLayout.col("umap_0").asArray() as NumberArray)[0],
      quantile([0.5], obsLayout.col("umap_1").asArray() as NumberArray)[0],
    ];

    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    centroidResult.forEach((coordinate: any) => {
      expect(coordinate).toEqual(expectedResult);
    });
  });

  test("field3 (boolean obsAnnotation)", () => {
    const centroidResult = calcCentroid(
      schema,
      "field3",
      obsAnnotations,
      { current: "umap", currentDimNames: ["umap_0", "umap_1"] },
      obsLayout
    );

    // Check to see that a centroid has been calculated for every categorical value
    const keysAsArray = Array.from(centroidResult.keys());
    expect(keysAsArray).toEqual(expect.arrayContaining([false, true]));

    // This expected result assumes that all cells belong in all categorical values inside of sample response
    const expectedResult = [
      quantile([0.5], obsLayout.col("umap_0").asArray() as NumberArray)[0],
      quantile([0.5], obsLayout.col("umap_1").asArray() as NumberArray)[0],
    ];

    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    centroidResult.forEach((coordinate: any) => {
      expect(coordinate).toEqual(expectedResult);
    });
  });
});
