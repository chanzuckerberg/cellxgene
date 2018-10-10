import _ from "lodash";
import * as Universe from "../../../src/util/stateManager/universe";
import * as World from "../../../src/util/stateManager/world";
import Crossfilter from "../../../src/util/typedCrossfilter";
import * as REST from "./sampleResponses";

/*
TODO: endpoints to test:

createWorldFromEntireUniverse(universe)
createWorldFromCurrentSelection(universe, world, crossfilter)
createObsDimensionMap(crossfilter, world)
subsetVarData(world, universe, varData)

*/

describe("createWorldFromEntireUniverse", () => {
  test("create from REST sample", () => {
    const universe = Universe.createUniverseFromRestV02Response(
      REST.config,
      REST.schema,
      REST.annotationsObs,
      REST.annotationsVar,
      REST.layoutObs
    );
    expect(universe).toBeDefined();

    const world = World.createWorldFromEntireUniverse(universe);
    expect(world).toBeDefined();

    expect(world).toMatchObject(
      expect.objectContaining({
        api: "0.2",
        nObs: universe.nObs,
        nVar: universe.nVar,
        schema: universe.schema,
        obsAnnotations: universe.obsAnnotations,
        varAnnotations: universe.varAnnotations,
        obsLayout: universe.obsLayout,

        summary: expect.objectContaining({
          obs: _(REST.schema.schema.annotations.obs)
            .keyBy("name")
            .mapValues(() => expect.any(Object))
            .value(),
          var: {} // TODO: currently unimplemneted
        }),

        varDataCache: expect.any(Object),

        worldObsIndex: null // indicating full universe
      })
    );
  });
});

describe("createWorldFromCurrentSelection", () => {
  test("create from REST sample", () => {
    /* create universe */
    const universe = Universe.createUniverseFromRestV02Response(
      REST.config,
      REST.schema,
      REST.annotationsObs,
      REST.annotationsVar,
      REST.layoutObs
    );
    /* create world */
    const originalWorld = World.createWorldFromEntireUniverse(universe);
    /* create crossfilter */
    const crossfilter = Crossfilter(originalWorld.obsAnnotations);

    /* fake a selection */
    const dimensionMap = World.createObsDimensionMap(
      crossfilter,
      originalWorld
    );
    dimensionMap.field1.filterRange(0, 3);
    dimensionMap.field3.filterExact(0);

    /* create the world from the selection */
    const world = World.createWorldFromCurrentSelection(
      universe,
      originalWorld,
      crossfilter
    );
    expect(world).toBeDefined();
  });
});
