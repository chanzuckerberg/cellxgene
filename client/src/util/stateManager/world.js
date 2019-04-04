// jshint esversion: 6

import { layoutDimensionName, obsAnnoDimensionName } from "../nameCreators";
import * as Dataframe from "../dataframe";

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
  return {
    /* schema/version related */
    schema: null,
    nObs: 0,
    nVar: 0,
    continuousPercentileMin: 0,
    continuousPercentileMax: 1,

    /* annotations */
    obsAnnotations: Dataframe.Dataframe.empty(),
    varAnnotations: Dataframe.Dataframe.empty(),

    /* layout of graph. Dataframe. */
    obsLayout: Dataframe.Dataframe.empty(),

    /*
    Var data columns - subset of all data (may be empty)
    */
    varData: Dataframe.Dataframe.empty(null, new Dataframe.KeyIndex())
  };
}

export function createWorldFromEntireUniverse(
  universe,
  continuousPercentileMin,
  continuousPercentileMax
) {
  const world = templateWorld();

  /*
  public interface follows
  */

  /* Schema related */
  world.schema = universe.schema;
  world.nObs = universe.nObs;
  world.nVar = universe.nVar;
  world.continuousPercentileMin = continuousPercentileMin;
  world.continuousPercentileMax = continuousPercentileMax;

  /* annotation dataframes */
  world.obsAnnotations = universe.obsAnnotations.clone();
  world.varAnnotations = universe.varAnnotations.clone();

  /* layout and display characteristics dataframe */
  world.obsLayout = universe.obsLayout.clone();

  /*
  Var data columns - subset of all
  */
  world.varData = universe.varData.clone();

  return world;
}

export function createWorldFromCurrentSelection(universe, world, crossfilter) {
export function createWorldFromCurrentSelection(
  universe,
  world,
  crossfilter,
  continuousPercentileMin,
  continuousPercentileMax
) {
  const newWorld = templateWorld();

  /* these don't change as only OBS are selected in our current implementation */
  newWorld.nVar = universe.nVar;
  newWorld.schema = universe.schema;
  newWorld.varAnnotations = universe.varAnnotations;

  /* now subset/cut obs */
  const mask = crossfilter.allSelectedMask();
  newWorld.obsAnnotations = world.obsAnnotations.isubsetMask(mask);
  newWorld.obsLayout = world.obsLayout.isubsetMask(mask);
  newWorld.nObs = newWorld.obsAnnotations.dims[0];

  newWorld.continuousPercentileMin = continuousPercentileMin;
  newWorld.continuousPercentileMax = continuousPercentileMax;

  /*
  Var data columns - subset of all
  */
  if (world.varData.isEmpty()) {
    newWorld.varData = world.varData.clone();
  } else {
    newWorld.varData = world.varData.isubsetMask(mask);
  }
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
