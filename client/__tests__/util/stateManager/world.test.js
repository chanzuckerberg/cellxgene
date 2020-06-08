import _ from "lodash";
import * as Universe from "../../../src/util/stateManager/universe";
import { matrixFBSToDataframe } from "../../../src/util/stateManager/matrix";
import * as World from "../../../src/util/stateManager/world";
import * as Dataframe from "../../../src/util/dataframe";
import Crossfilter from "../../../src/util/typedCrossfilter";
import { DimTypes } from "../../../src/util/typedCrossfilter/crossfilter";
import * as REST from "./sampleResponses";
import {
  obsAnnoDimensionName,
  layoutDimensionName,
} from "../../../src/util/nameCreators";

/*
Helper - creates universe, world, corssfilter and dimensionMap from
the default REST test response.
*/
const defaultBigBang = () => {
  /* create unverse, world, crossfilter and dimensionMap */
  /* create universe */
  let universe = Universe.createUniverseFromResponse(
    _.cloneDeep(REST.config),
    _.cloneDeep(REST.schema)
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

  /* create world */
  const world = World.createWorldFromEntireUniverse(universe);
  /* create crossfilter */
  const crossfilter = World.createObsDimensions(
    new Crossfilter(world.obsAnnotations),
    world,
    REST.schema.schema.layout.obs[0].dims
  );

  return {
    universe,
    world,
    crossfilter,
  };
};

describe("createWorldFromEntireUniverse", () => {
  test("create from REST sample", () => {
    const universe = Universe.createUniverseFromResponse(
      _.cloneDeep(REST.config),
      _.cloneDeep(REST.schema),
      matrixFBSToDataframe(_.cloneDeep(REST.annotationsObs)),
      matrixFBSToDataframe(_.cloneDeep(REST.annotationsVar)),
      matrixFBSToDataframe(_.cloneDeep(REST.layoutObs))
    );
    expect(universe).toBeDefined();

    const world = World.createWorldFromEntireUniverse(universe);
    expect(world).toBeDefined();

    expect(world).toMatchObject(
      expect.objectContaining({
        nObs: universe.nObs,
        nVar: universe.nVar,
        schema: universe.schema,
        obsAnnotations: expect.any(Dataframe.Dataframe),
        varAnnotations: expect.any(Dataframe.Dataframe),
        obsLayout: expect.any(Dataframe.Dataframe),
        varData: expect.any(Dataframe.Dataframe),
        clipQuantiles: { min: 0, max: 1 },
        unclipped: {
          obsAnnotations: expect.any(Dataframe.Dataframe),
          varData: expect.any(Dataframe.Dataframe),
        },
      })
    );
  });
});

describe("createWorldFromCurrentSelection", () => {
  test("create from REST sample", () => {
    const {
      universe,
      world: originalWorld,
      crossfilter: originalCrossfilter,
    } = defaultBigBang();

    /* mock a selection */
    const crossfilter = originalCrossfilter
      .select(obsAnnoDimensionName("field1"), { mode: "range", lo: 0, hi: 5 })
      .select(obsAnnoDimensionName("field3"), {
        mode: "exact",
        values: [false],
      });

    /* create the world from the selection */
    const world = World.createWorldBySelection(
      universe,
      originalWorld,
      crossfilter
    );
    expect(world).toBeDefined();
    expect(world.nObs).toEqual(crossfilter.countSelected());

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
      .filter((idx) => matchFilter(universe.obsAnnotations, idx))
      .value();

    expect(world).toMatchObject(
      expect.objectContaining({
        nObs: matchingIndices.length,
        nVar: universe.nVar,
        schema: universe.schema,
        clipQuantiles: { min: 0, max: 1 },
        obsAnnotations: expect.any(Dataframe.Dataframe),
        varAnnotations: expect.any(Dataframe.Dataframe),
        obsLayout: expect.any(Dataframe.Dataframe),
        varData: expect.any(Dataframe.Dataframe),
        unclipped: {
          obsAnnotations: expect.any(Dataframe.Dataframe),
          varData: expect.any(Dataframe.Dataframe),
        },
      })
    );

    expect(world.obsAnnotations.rowIndex.labels()).toEqual(
      new Int32Array(matchingIndices)
    );
    expect(world.obsAnnotations.colIndex.labels()).toEqual(
      universe.obsAnnotations.colIndex.labels()
    );
    expect(world.obsLayout.rowIndex.labels()).toEqual(
      new Int32Array(matchingIndices)
    );
    expect(world.obsLayout.colIndex.labels()).toEqual(
      world.schema.layout.obs[0].dims
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

    const { crossfilter } = defaultBigBang();
    const annotationNames = _.map(
      REST.schema.schema.annotations.obs.columns,
      (c) => c.name
    );
    const obsIndexColName = REST.schema.schema.annotations.obs.index;
    const schemaByObsName = _.keyBy(
      REST.schema.schema.annotations.obs.columns,
      "name"
    );
    expect(crossfilter).toBeDefined();
    annotationNames.forEach((name) => {
      const dim = crossfilter.dimensions[obsAnnoDimensionName(name)];
      if (name === obsIndexColName) {
        expect(dim).toBeUndefined();
      } else {
        const { type } = schemaByObsName[name];
        if (type === "string" || type === "boolean" || type === "categorical") {
          expect(dim.dim).toBeInstanceOf(DimTypes.enum);
        } else {
          expect(dim.dim).toBeInstanceOf(DimTypes.scalar);
        }
      }
    });
    expect(
      crossfilter.dimensions[layoutDimensionName("XY")].dim
    ).toBeInstanceOf(DimTypes.spatial);
  });
});

describe("worldEqUniverse", () => {
  const { universe, world } = defaultBigBang();
  const result = World.worldEqUniverse(world, universe);
  expect(result).toBe(true);
});
