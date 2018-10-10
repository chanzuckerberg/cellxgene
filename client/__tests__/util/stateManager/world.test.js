import _ from "lodash";
import * as Universe from "../../../src/util/stateManager/universe";
import * as World from "../../../src/util/stateManager/world";
import Crossfilter from "../../../src/util/typedCrossfilter";
import * as REST from "./sampleResponses";

/*
TODO: endpoints to test:

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
    dimensionMap.field1.filterRange([0, 5]);
    dimensionMap.field3.filterExact(false);

    /* create the world from the selection */
    const world = World.createWorldFromCurrentSelection(
      universe,
      originalWorld,
      crossfilter
    );
    expect(world).toBeDefined();
    expect(world.nObs).toEqual(crossfilter.countFiltered());

    /*
    calculate expected values and match against result
    */

    /* matchFilter must match the dimension filters above */
    const matchFilter = val => val.field1 >= 0 && val.field1 < 5 && !val.field3;
    const universeIndices = _()
      .range(universe.nObs)
      .filter(idx => matchFilter(universe.obsAnnotations[idx]))
      .value();

    const expected = {
      nObs: universeIndices.length,
      obsAnnotations: _.map(universeIndices, i => universe.obsAnnotations[i]),
      obsLayout: {
        X: _.map(universeIndices, i => universe.obsLayout.X[i]),
        Y: _.map(universeIndices, i => universe.obsLayout.Y[i])
      },
      worldObsIndex: _.transform(
        universeIndices,
        (result, univIdx, worldIdx) => {
          result[univIdx] = worldIdx;
        },
        new Array(universe.nObs).fill(-1)
      )
    };

    expect(world).toMatchObject(
      expect.objectContaining({
        api: "0.2",
        nObs: expected.nObs,
        nVar: universe.nVar,
        schema: universe.schema,
        obsAnnotations: expected.obsAnnotations,
        varAnnotations: universe.varAnnotations,
        obsLayout: expected.obsLayout,
        summary: expect.any(Object) /* we could do better! */,
        varDataCache: expect.any(Object),
        worldObsIndex: expected.worldObsIndex
      })
    );
  });
});
