// jshint esversion: 6

import _ from "lodash";
import * as kvCache from "./keyvalcache";
import summarizeAnnotations from "./summarizeAnnotations";
import { layoutDimensionName, obsAnnoDimensionName } from "../nameCreators";
import Crossfilter from "../typedCrossfilter";
import { sliceByIndex } from "../typedCrossfilter/util";

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

* summary: summary of each obsAnnotation column (eg, numeric extent for
  continuous data, category counts for categorical metadata)

* varDataCache: expression columns, in a kvCache.   TODO: maybe move to a
  Dataframe in the future.

*/

/* varDataCache config - see kvCache for semantics */
const VarDataCacheLowWatermark = 32; // cache element count
const VarDataCacheTTLMs = 1000; // min cache time in MS

function templateWorld() {
  return {
    /* schema/version related */
    api: null,
    schema: null,
    nObs: 0,
    nVar: 0,

    /* annotations */
    obsAnnotations: null,
    varAnnotations: null,

    /* layout of graph. Dataframe. */
    obsLayout: null,

    /* derived data summaries XXX: consider exploding in place */
    summary: null,

    varDataCache: kvCache.create(
      VarDataCacheLowWatermark,
      VarDataCacheTTLMs
    ) /* cache of var data (expression) */
  };
}

export function createWorldFromEntireUniverse(universe) {
  if (!universe.finalized) {
    throw new Error("World can't be created from an partial Universe");
  }

  const world = templateWorld();

  /*
  public interface follows
  */

  /* Schema related */
  world.api = universe.api;
  world.schema = universe.schema;
  world.nObs = universe.nObs;
  world.nVar = universe.nVar;

  /* annotations */
  world.obsAnnotations = universe.obsAnnotations;
  world.varAnnotations = universe.varAnnotations;

  /* layout and display characteristics */
  world.obsLayout = universe.obsLayout;

  /* derived data & summaries */
  world.summary = summarizeAnnotations(
    world.schema,
    world.obsAnnotations,
    world.varAnnotations
  );

  /* build the varDataCache */
  world.varDataCache = kvCache.map(
    universe.varDataCache,
    val => subsetVarData(world, universe, val),
    { lowWatermark: VarDataCacheLowWatermark, minTTL: VarDataCacheTTLMs }
  );

  return world;
}

export function createWorldFromCurrentSelection(universe, world, crossfilter) {
  const newWorld = templateWorld();

  /* these don't change as only OBS are selected in our current implementation */
  newWorld.api = universe.api;
  newWorld.nVar = universe.nVar;
  newWorld.schema = universe.schema;
  newWorld.varAnnotations = universe.varAnnotations;

  /* now subset/cut obs */
  const mask = crossfilter.allFilteredMask();
  newWorld.obsAnnotations = world.obsAnnotations.icutByMask(mask);
  newWorld.obsLayout = world.obsLayout.icutByMask(mask);
  newWorld.nObs = newWorld.obsAnnotations.dims[0];

  /* derived data & summaries */
  newWorld.summary = summarizeAnnotations(
    newWorld.schema,
    newWorld.obsAnnotations,
    newWorld.varAnnotations
  );

  /* build the varDataCache */
  newWorld.varDataCache = kvCache.map(
    universe.varDataCache,
    val => subsetVarData(newWorld, universe, val),
    { lowWatermark: VarDataCacheLowWatermark, minTTL: VarDataCacheTTLMs }
  );
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

/*
  Return a crossfilter dimension for the specified world & named gene.

  NOTE: this assumes that the expression data was already loaded,
  by calling an appropriate action creator.

  Caller needs to *save* this dimension somewhere for it to be later used.
  Dimension must be destroyed by calling dimension.dispose()
  when it is no longer needed
  (it will not be garbage collected without this call)
*/

export function createVarDimension(
  world,
  _worldVarDataCache,
  crossfilter,
  geneName
) {
  // return crossfilter.dimension(_worldVarDataCache[geneName], Float32Array);
  return crossfilter.dimension(
    Crossfilter.ScalarDimension,
    _worldVarDataCache[geneName],
    Float32Array
  );
}

export function createObsDimensionMap(crossfilter, world) {
  /*
  create and return a crossfilter dimension for every obs annotation
  for which we have a supported type.
  */
  const { schema, obsLayout, obsAnnotations } = world;

  // Create a crossfilter dimension for all obs annotations *except* 'name'
  const dimensionMap = _(schema.annotations.obs)
    .filter(anno => anno.name !== "name")
    .transform((result, anno) => {
      const dimType = deduceDimensionType(anno, anno.name);
      const colData = obsAnnotations.col(anno.name).asArray();
      if (dimType === "enum") {
        result[obsAnnoDimensionName(anno.name)] = crossfilter.dimension(
          Crossfilter.EnumDimension,
          colData
        );
      } else if (dimType) {
        result[obsAnnoDimensionName(anno.name)] = crossfilter.dimension(
          Crossfilter.ScalarDimension,
          colData,
          dimType
        );
      } // else ignore the annotation
    }, {})
    .value();

  /*
  Add crossfilter dimensions allowing filtering on layout
  */
  dimensionMap[layoutDimensionName("XY")] = crossfilter.dimension(
    Crossfilter.SpatialDimension,
    obsLayout.col("X").asArray(),
    obsLayout.col("Y").asArray()
  );

  return dimensionMap;
}

export function worldEqUniverse(world, universe) {
  return world.obsAnnotations === universe.obsAnnotations;
}

export function subsetVarData(world, universe, varData) {
  // If world === universe, just return the entire varData array
  if (worldEqUniverse(world, universe)) {
    return varData;
  }
  return sliceByIndex(varData, world.obsAnnotations.rowIndex.keys());
}

export function getSelectedByIndex(crossfilter) {
  /*
  return array of obsIndex, containing all selected obs/cells.
  */
  const selected = crossfilter.allFilteredMask(); // array of bool-ish
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
