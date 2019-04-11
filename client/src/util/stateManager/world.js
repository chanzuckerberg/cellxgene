// jshint esversion: 6

import { layoutDimensionName, obsAnnoDimensionName } from "../nameCreators";
import * as Dataframe from "../dataframe";
import ImmutableTypedCrossfilter from "../typedCrossfilter/crossfilter";

/*

World is a subset of universe.   Most code should use world, and should
(generally) not use Universe.   World contains any per-obs or per-var data
that must be consistent acorss the app when we view/manipulate subsets
of Universe.

Private API indicated by leading underscore in key name (eg, _foo).  Anything else
is public.

Notable keys in the world object:

* nObs, nVar: dimensions

* schema: data schema from the server

* clipQuantiles: the quantiles used to clip all data in world.

* obsAnnotations:

  Dataframe containing obs annotations.  Columns are indexed by annotation
  name (eg, 'tissue type'), and rows are indexed by the REST API obsIndex
  (ie, the offset into the underlying server-side dataframe).

  This indexing means that you can access data by _either_ the server's
  obxIndex, or the offset into the client-side column array .  Be careful
  to know which you want and are using.

* obsLayout:

  A dataframe containing the X/Y layout for all obs.  Columns are named
  'X' and 'Y', and rows are indexed in the same way as obsAnnotation.

* varData: a cache of expression columns, stored in a Dataframe.  Cache
  managed by controls reducer.

*/

function templateWorld() {
  const obsAnnotations = Dataframe.Dataframe.empty();
  const varData = Dataframe.Dataframe.empty(null, new Dataframe.KeyIndex());
  return {
    /* schema/version related */
    schema: null,
    nObs: 0,
    nVar: 0,
    clipQuantiles: { min: 0, max: 1 },

    /* annotations */
    obsAnnotations,
    varAnnotations: Dataframe.Dataframe.empty(),

    /* layout of graph. Dataframe. */
    obsLayout: Dataframe.Dataframe.empty(),

    /* Var data columns - subset of all data (may be empty) */
    varData,

    /* unclipped dataframes - subset, but not value clipped */
    unclipped: {
      obsAnnotations,
      varData
    }
  };
}

function clipDataframe(df, lowerQuantile, upperQuantile, value = Number.NaN) {
  /*
  Clip all values above or below specified quantiles to `value`.
  */

  // XXX - this is not finished

  if (lowerQuantile < 0) lowerQuantile = 0;
  if (upperQuantile > 1) upperQuantile = 1;
  if (lowerQuantile === 0 && upperQuantile === 1) return df;

  const keys = df.keys();
  return df.mapColumns((col, colIdx) => {
    /* do we actually have access to the colIdx from `mapColumns`?*/
    const colName = keys[colIdx];
    /* how do we check if it's continuous metadata?
    if (!colName_is_a_continuous_metadata_field) return col;
    */

    const newCol = col.slice();
    for (let i = 0, l = newCol.length; i < l; i += 1) {
      const colMin = ImmutableTypedCrossfilter.percentile(
        colName,
        lowerQuantile
      );
      const colMax = ImmutableTypedCrossfilter.percentile(
        colName,
        upperQuantile
      );
      if (newCol[i] < colMin || newCol[i] > colMax) {
        newCol[i] = value;
      }
    }
    return newCol;
  });
}

/*
Create World with contents eq entire universe.   Commonly used to initialize World.
If clipQuantiles
*/
export function createWorldFromEntireUniverse(universe) {
  const world = templateWorld();

  /* Schema related */
  world.schema = universe.schema;
  world.nObs = universe.nObs;
  world.nVar = universe.nVar;
  world.clipQuantiles = { min: 0, max: 1 };

  /* dataframes: annotations and layout */
  world.obsAnnotations = universe.obsAnnotations.clone();
  world.varAnnotations = universe.varAnnotations.clone();
  world.obsLayout = universe.obsLayout.clone();

  /* Var dataframe - contains a subset of all var columns */
  world.varData = universe.varData.clone();

  world.unclipped = {
    obsAnnotations: world.obsAnnotations,
    varData: world.varData
  };

  return world;
}

/*
Create world as a subset of the current world.   Params:
  * universe: the universe object
  * world: the _current_ world object
  * crossfilter: If crossfilter is specified, the new world will contain
    only those elements currently selected.
  * clipQuantiles: if specified, the new worlds values will be clipped.

TODO/XXX: this function needs a better name.
*/
export function createWorldFromCurrentWorld(
  universe,
  world,
  crossfilter = null,
  clipQuantiles
) {
  const newWorld = templateWorld();

  /* these don't change as only OBS are selected in our current implementation */
  newWorld.nVar = universe.nVar;
  newWorld.schema = universe.schema;
  newWorld.varAnnotations = universe.varAnnotations;

  /* if not specified, inherit clip quantiles from current World */
  if (!clipQuantiles) {
    newWorld.clipQuantiles = world.clipQuantiles;
  } else {
    newWorld.clipQuantiles = clipQuantiles;
  }

  /* if crossfilter provided, subset/cut world to current selection */
  if (!crossfilter) {
    newWorld.obsAnnotations = world.obsAnnotations.clone();
    newWorld.obsLayout = world.obsLayout.clone();
    newWorld.varData = world.varData.clone();
  } else {
    const mask = crossfilter.allSelectedMask();
    newWorld.obsAnnotations = world.obsAnnotations.isubsetMask(mask);
    newWorld.obsLayout = world.obsLayout.isubsetMask(mask);

    /* Var data columns - subset of all */
    if (world.varData.isEmpty()) {
      newWorld.varData = world.varData.clone();
    } else {
      newWorld.varData = world.varData.isubsetMask(mask);
    }
  }

  newWorld.nObs = newWorld.obsAnnotations.dims[0];
  return newWorld;
}

/*
 Deduce the correct crossfilter dimension type from a metadata
 schema description.
*/
function deduceDimensionType(attributes, fieldName) {
  let dimensionType;
  const { type } = attributes;
  if (type === "string" || type === "categorical" || type === "boolean") {
    dimensionType = "enum";
  } else if (type === "int32") {
    dimensionType = Int32Array;
  } else if (type === "float32") {
    dimensionType = Float32Array;
  } else {
    /*
    Currently not supporting boolean and categorical types.
    */
    console.error(
      `Warning - REST API returned unknown metadata schema (${type}) for field ${fieldName}.`
    );
    // skip it - we don't know what to do with this type
  }
  return dimensionType;
}

export function createObsDimensions(crossfilter, world) {
  /*
  create and return a crossfilter with a dimension for every obs annotation
  for which we have a supported type, *except* 'name'
  */
  const { schema, obsLayout, obsAnnotations } = world;
  const annoList = schema.annotations.obs.filter(anno => anno.name !== "name");
  crossfilter = annoList.reduce((xfltr, anno) => {
    const dimType = deduceDimensionType(anno, anno.name);
    const colData = obsAnnotations.col(anno.name).asArray();
    const name = obsAnnoDimensionName(anno.name);
    if (dimType === "enum") {
      return xfltr.addDimension(name, "enum", colData);
    }
    if (dimType) {
      return xfltr.addDimension(name, "scalar", colData, dimType);
    }
    return xfltr;
  }, crossfilter);

  return crossfilter.addDimension(
    layoutDimensionName("XY"),
    "spatial",
    obsLayout.col("X").asArray(),
    obsLayout.col("Y").asArray()
  );
}

export function worldEqUniverse(world, universe) {
  return (
    world.obsAnnotations === universe.obsAnnotations ||
    world.obsAnnotations.rowIndex === universe.obsAnnotations.rowIndex
  );
}

export function getSelectedByIndex(crossfilter) {
  /*
  return array of obsIndex, containing all selected obs/cells.
  */
  const selected = crossfilter.allSelectedMask(); // array of bool-ish
  const keys = crossfilter.data.rowIndex.keys(); // row keys, aka universe rowIndex

  const set = new Int32Array(selected.length);
  let numElems = 0;
  for (let i = 0, l = selected.length; i < l; i += 1) {
    if (selected[i]) {
      set[numElems] = keys[i];
      numElems += 1;
    }
  }
  return new Int32Array(set.buffer, 0, numElems);
}
