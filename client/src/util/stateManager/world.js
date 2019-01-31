// jshint esversion: 6

import _ from "lodash";
import * as kvCache from "./keyvalcache";
import summarizeAnnotations from "./summarizeAnnotations";
import { layoutDimensionName, obsAnnoDimensionName } from "../nameCreators";
import { sliceByIndex } from "../typedCrossfilter/util";

/*  XXX comments need revising to accomadate dataframe changes

World is a subset of universe.   Most code should use world, and should
(generally) not use Universe.   World contains any per-obs or per-var data
that must be consistent acorss the app when we view/manipulate subsets
of Universe.

Private API indicated by leading underscore in key name (eg, _foo).  Anything else
is public.

World contains several public keys, obsAnnotations, and obsLayout, which are
dataframes contianing information about an OBS in the same order/offset.

* obsAnnotations:

  obsAnnotations will return an array of objects.  Each object contains all annotation
  values for a given observation/cell, keyed by annotation name, PLUS a key
  '__cellId__', containing a REST API ID for this obs/cell (referred to as the
  obsIndex in the REST 0.2 spec or cellIndex in the 0.1 spec.

  Example:  [ { __cellId__: 99, cluster: 'blue', numReads: 93933 } ]

  NOTE: world.obsAnnotation should be identical to the old state.cells value,
  EXCEPT that
    * __cellIndex__ renamed to __index__
    * __x__ and __y__ are now in world.obsLayout
    * __color__ and __colorRBG__ should be moved to controls reducer

* obsLayout:

  obsLayout will return an object containing two arrays, containing X and Y
  coordinates respectively.

  Example: { X: [ 0.33, 0.23, ... ], Y: [ 0.8, 0.777, ... ]}

* crossfilter - a crossfilter object across world.obsAnnotations

* dimensionMap - an object mapping annotation names to dimensions on
  the crossfilter

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
  newWorld.obsAnnotations = universe.obsAnnotations.icutByMask(mask, null);
  newWorld.obsLayout = universe.obsLayout.icutByMask(mask, null);
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
  return crossfilter.dimension(_worldVarDataCache[geneName], Float32Array);
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
      if (dimType) {
        const colData = obsAnnotations.col(anno.name).asArray();
        result[obsAnnoDimensionName(anno.name)] = crossfilter.dimension(
          colData,
          dimType
        );
      } // else ignore the annotation
    }, {})
    .value();

  /*
  Add crossfilter dimensions allowing filtering on layout
  */
  const X = obsLayout.col("X").asArray();
  const Y = obsLayout.col("Y").asArray();
  dimensionMap[layoutDimensionName("X")] = crossfilter.dimension(
    X,
    X.constructor
  );
  dimensionMap[layoutDimensionName("Y")] = crossfilter.dimension(
    Y,
    Y.constructor
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
