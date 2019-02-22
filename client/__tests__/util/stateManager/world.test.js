import _ from "lodash";
import * as Universe from "../../../src/util/stateManager/universe";
import * as World from "../../../src/util/stateManager/world";
import * as Dataframe from "../../../src/util/dataframe";
import Crossfilter from "../../../src/util/typedCrossfilter";
import * as REST from "./sampleResponses";
import {
  obsAnnoDimensionName,
  layoutDimensionName
} from "../../../src/util/nameCreators";

/*
Helper - creates universe, world, corssfilter and dimensionMap from
the default REST test response.
*/
const defaultBigBang = () => {
  /* create unverse, world, crossfilter and dimensionMap */
  /* create universe */
  const universe = Universe.createUniverseFromResponse(
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
    const universe = Universe.createUniverseFromResponse(
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
        nObs: universe.nObs,
        nVar: universe.nVar,
        schema: universe.schema,
        obsAnnotations: universe.obsAnnotations,
        varAnnotations: universe.varAnnotations,
        obsLayout: universe.obsLayout,
        varData: expect.any(Dataframe.Dataframe)
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
    const matchFilter = (df, row) => {
      const field1 = df.at(row, "field1");
      const field3 = df.at(row, "field3");
      return field1 >= 0 && field1 < 5 && !field3;
    };
    const matchingIndices = _()
      .range(universe.nObs)
      .filter(idx => matchFilter(universe.obsAnnotations, idx))
      .value();

    expect(world).toMatchObject(
      expect.objectContaining({
        nObs: matchingIndices.length,
        nVar: universe.nVar,
        schema: universe.schema,
        obsAnnotations: expect.any(Dataframe.Dataframe),
        varAnnotations: universe.varAnnotations,
        obsLayout: expect.any(Dataframe.Dataframe),
        varData: expect.any(Dataframe.Dataframe)
      })
    );

    expect(world.obsAnnotations.rowIndex.keys()).toEqual(
      new Int32Array(matchingIndices)
    );
    expect(world.obsAnnotations.colIndex.keys()).toEqual(
      universe.obsAnnotations.colIndex.keys()
    );
    expect(world.obsLayout.rowIndex.keys()).toEqual(
      new Int32Array(matchingIndices)
    );
    expect(world.obsLayout.colIndex.keys()).toEqual(["X", "Y"]);
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
    expect(newWorld.obsAnnotations.rowIndex.keys()).toEqual(
      new Int32Array([0, 2])
    );

    /* expect a subset */
    const result = World.subsetVarData(newWorld, universe, sourceVarData);
    expect(result).not.toBe(sourceVarData);
    expect(result).toHaveLength(newWorld.nObs);
    /* check that we have expected source var content */
    expect(result).toMatchObject(new Float32Array([0, 2]));
  });
});

describe("createVarDataDimension", () => {
  /* create default universe */
  const { world, crossfilter } = defaultBigBang();
  world.varData = world.varData.withCol(
    "GENE",
    Float32Array.from(_.range(world.nObs))
  );
  const result = World.createVarDataDimension(world, crossfilter, "GENE");
  expect(result).toBeInstanceOf(Crossfilter.ScalarDimension);
});

describe("worldEqUniverse", () => {
  const { universe, world } = defaultBigBang();
  const result = World.worldEqUniverse(world, universe);
  expect(result).toBe(true);
});
