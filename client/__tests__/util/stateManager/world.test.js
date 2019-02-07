import _ from "lodash";
import * as Universe from "../../../src/util/stateManager/universe";
import * as World from "../../../src/util/stateManager/world";
import Crossfilter from "../../../src/util/typedCrossfilter";
import * as REST from "./sampleResponses";
import {
  obsAnnoDimensionName,
  layoutDimensionName
} from "../../../src/util/nameCreators";
import * as kvCache from "../../../src/util/stateManager/keyvalcache";

/*
Helper - creates universe, world, corssfilter and dimensionMap from
the default REST test response.
*/
const defaultBigBang = () => {
  /* create unverse, world, crossfilter and dimensionMap */
  /* create universe */
  const universe = Universe.createUniverseFromRestV02Response(
    REST.config,
    REST.schema,
    REST.annotationsObs,
    REST.annotationsVar,
    REST.layoutObs
  );
  /* create world */
  const world = World.createWorldFromEntireUniverse(universe);
  /* create crossfilter */
  const crossfilter = Crossfilter(world.obsAnnotations);
  /* create dimension map */
  const dimensionMap = World.createObsDimensionMap(crossfilter, world);

  return {
    universe,
    world,
    crossfilter,
    dimensionMap
  };
};

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
            .filter(v => v.name !== "name")
            .keyBy("name")
            .mapValues(() => expect.any(Object))
            .value(),
          var: _(REST.schema.schema.annotations.var)
            .filter(v => v.name !== "name")
            .keyBy("name")
            .mapValues(() => expect.any(Object))
            .value()
        }),

        varDataCache: expect.any(Object),

        obsIndex: null, // null indicating full universe
        obsBackIndex: null
      })
    );
  });
});

describe("createWorldFromCurrentSelection", () => {
  test("create from REST sample", () => {
    const {
      universe,
      world: originalWorld,
      crossfilter,
      dimensionMap
    } = defaultBigBang();

    /* mock a selection */
    dimensionMap[obsAnnoDimensionName("field1")].filterRange([0, 5]);
    dimensionMap[obsAnnoDimensionName("field3")].filterExact(false);

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
        X: new Float32Array(
          _.map(universeIndices, i => universe.obsLayout.X[i])
        ),
        Y: new Float32Array(
          _.map(universeIndices, i => universe.obsLayout.Y[i])
        )
      },
      obsBackIndex: _.transform(
        universeIndices,
        (result, univIdx, worldIdx) => {
          result[univIdx] = worldIdx;
        },
        new Uint32Array(universe.nObs).fill(-1)
      ),
      obsIndex: new Uint32Array(universeIndices)
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
        summary: {
          obs: expect.any(Object) /* we could do better! */,
          var: expect.any(Object) /* we could do better! */
        },
        varDataCache: expect.any(Object),
        obsIndex: expected.obsIndex,
        obsBackIndex: expected.obsBackIndex
      })
    );
  });
});

describe("createObsDimensionMap", () => {
  test("when universe eq world", () => {
    /*
    check for:
    - creates a dimension for all obsAnnotations, PLUS X/Y layout
    - check that dimension typing is sane
    */

    const { dimensionMap } = defaultBigBang();
    const annotationNames = _.map(
      REST.schema.schema.annotations.obs,
      c => c.name
    );
    const schemaByObsName = _.keyBy(REST.schema.schema.annotations.obs, "name");
    expect(dimensionMap).toBeDefined();
    annotationNames.forEach(name => {
      const dim = dimensionMap[obsAnnoDimensionName(name)];
      if (name === "name") {
        expect(dim).toBeUndefined();
      } else {
        const { type } = schemaByObsName[name];
        if (type === "string" || type === "boolean" || type === "categorical") {
          expect(dim).toBeInstanceOf(Crossfilter.EnumDimension);
        } else {
          expect(dim).toBeInstanceOf(Crossfilter.ScalarDimension);
        }
      }
    });
    expect(dimensionMap[layoutDimensionName("XY")]).toBeInstanceOf(
      Crossfilter.SpatialDimension
    );
  });
});

describe("subsetVarData", () => {
  test("when world eq universe", () => {
    const { universe, world } = defaultBigBang();
    /* create a mock varData array for subsetting */
    const sourceVarData = new Float32Array(universe.nObs);

    /* expect literally the same object back */
    const result = World.subsetVarData(world, universe, sourceVarData);
    expect(result).toBe(sourceVarData);
  });

  test("when world neq universe", () => {
    const { universe, world, crossfilter, dimensionMap } = defaultBigBang();
    /* create a mock varData array for subsetting */
    const sourceVarData = Float32Array.from(_.range(universe.nObs));

    /* mock a selection */
    dimensionMap[obsAnnoDimensionName("field1")].filterRange([0, 5]);
    dimensionMap[obsAnnoDimensionName("field3")].filterExact(false);

    /* create the world from the selection */
    const newWorld = World.createWorldFromCurrentSelection(
      universe,
      world,
      crossfilter
    );
    expect(newWorld.obsIndex).toMatchObject(new Uint32Array([0, 2]));

    /* expect a subset */
    const result = World.subsetVarData(newWorld, universe, sourceVarData);
    expect(result).not.toBe(sourceVarData);
    expect(result).toHaveLength(newWorld.nObs);
    /* check that we have expected source var content */
    expect(result).toMatchObject(new Float32Array([0, 2]));
  });
});

describe("createVarDimension", () => {
  /* create default universe */
  const { world, crossfilter } = defaultBigBang();
  /* create a mock var data cache */
  const varDataCache = kvCache.set(
    kvCache.create(),
    "GENE",
    Float32Array.from(_.range(world.nObs))
  );
  const result = World.createVarDimension(
    world,
    varDataCache,
    crossfilter,
    "GENE"
  );
  expect(result).toBeInstanceOf(Crossfilter.ScalarDimension);
});

describe("worldEqUniverse", () => {
  const { universe, world } = defaultBigBang();
  const result = World.worldEqUniverse(world, universe);
  expect(result).toBe(true);
});
