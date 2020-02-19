import _ from "lodash";

import calcCentroid from "../../src/util/centroid";

import * as Universe from "../../src/util/stateManager/universe";
import * as World from "../../src/util/stateManager/world";
import * as REST from "./stateManager/sampleResponses";
import { ControlsHelpers as CH } from "../../src/util/stateManager";

describe("centroid", () => {
  let world;
  let categoricalSelection;
  beforeAll(() => {
    let universe = Universe.createUniverseFromResponse(
      _.cloneDeep(REST.config),
      _.cloneDeep(REST.schema)
    );

    universe = {
      ...universe,
      ...Universe.addObsAnnotations(
        universe,
        Universe.matrixFBSToDataframe(REST.annotationsObs)
      ),
      ...Universe.addVarAnnotations(
        universe,
        Universe.matrixFBSToDataframe(REST.annotationsVar)
      ),
      ...Universe.addObsLayout(
        universe,
        Universe.matrixFBSToDataframe(REST.layoutObs)
      )
    };

    /* create world */
    world = World.createWorldFromEntireUniverse(universe);

    categoricalSelection = CH.createCategoricalSelection(
      world,
      CH.selectableCategoryNames(world.schema, CH.maxCategoryItems(REST.config))
    );
  });
  test("memoization", () => {
    console.log("world", world);
    console.log("categoricalSelection", categoricalSelection);
  });
});
