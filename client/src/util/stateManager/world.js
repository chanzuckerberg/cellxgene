// jshint esversion: 6

import _ from "lodash";
import * as kvCache from "./keyvalcache";
import summarizeAnnotations from "./summarizeAnnotations";
import { layoutDimensionName, obsAnnoDimensionName } from "../nameCreators";

/*
World is a subset of universe.   Most code should use world, and should
(generally) not use Universe.   World contains any per-obs or per-var data
that must be consistent acorss the app when we view/manipulate subsets
of Universe.

Private API indicated by leading underscore in key name (eg, _foo).  Anything else
is public.

World contains several public keys, obsAnnotations, and obsLayout, which are
arrays contianing information about an OBS in the same order/offset.  In
other words, world.obsAnnotations[0] and world.obsLayout.X[0] refer to the same
obs/cell.

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
    // map from universe obsIndex to world offset.
    // Undefined / null indicates identity mapping.
    worldObsIndex: null,

    /* schema/version related */
    api: null,
    schema: null,
    nObs: 0,
    nVar: 0,

    /* annotations */
    obsAnnotations: null,
    varAnnotations: null,

    /* layout of graph */
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

  // map from the universe obsIndex to our world offset.
  // undefined/null indicates identity map.
  world.worldObsIndex = null;

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
  newWorld.api = world.api;
  newWorld.nVar = world.nVar;
  newWorld.schema = world.schema;
  newWorld.varAnnotations = world.varAnnotations;

  /*
  Subset world from universe based upon world's current selection.  Only those
  fields which are subset by observation selection/filtering need to be updated.
  */
  const numSelected = crossfilter.countFiltered();

  /*
  Create a world which is based upon current selection
  */
  newWorld.nObs = numSelected;
  newWorld.obsAnnotations = new Array(numSelected);
  newWorld.obsLayout = {
    X: new Array(numSelected),
    Y: new Array(numSelected)
  };
  newWorld.worldObsIndex = new Array(universe.nObs);

  for (let i = 0, sel = 0; i < world.nObs; i += 1) {
    if (crossfilter.isElementFiltered(i)) {
      newWorld.obsAnnotations[sel] = world.obsAnnotations[i];
      newWorld.obsLayout.X[sel] = world.obsLayout.X[i];
      newWorld.obsLayout.Y[sel] = world.obsLayout.Y[i];
      sel += 1;
    }
  }

  // build index to our world offset
  newWorld.worldObsIndex.fill(-1); // default - aka unused
  for (let i = 0; i < newWorld.nObs; i += 1) {
    newWorld.worldObsIndex[newWorld.obsAnnotations[i].__index__] = i;
  }

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
  const { worldObsIndex } = world;
  const varData = _worldVarDataCache[geneName];
  const worldIndex = worldObsIndex ? idx => worldObsIndex[idx] : idx => idx;

  return crossfilter.dimension(
    r => varData[worldIndex(r.__index__)],
    Float32Array
  );
}

export function createObsDimensionMap(crossfilter, world) {
  /*
  create and return a crossfilter dimension for every obs annotation
  for which we have a supported type.
  */
  const { schema, obsLayout, worldObsIndex } = world;

  const dimensionMap = _.transform(
    schema.annotations.obs,
    (result, anno) => {
      const dimType = deduceDimensionType(anno, anno.name);
      if (dimType) {
        result[obsAnnoDimensionName(anno.name)] = crossfilter.dimension(
          r => r[anno.name],
          dimType
        );
      } // else ignore the annotation
    },
    {}
  );

  /*
  Add crossfilter dimensions allowing filtering on layout
  */
  const worldIndex = worldObsIndex ? idx => worldObsIndex[idx] : idx => idx;
  dimensionMap[layoutDimensionName("X")] = crossfilter.dimension(
    r => obsLayout.X[worldIndex(r.__index__)],
    Float32Array
  );
  dimensionMap[layoutDimensionName("Y")] = crossfilter.dimension(
    r => obsLayout.Y[worldIndex(r.__index__)],
    Float32Array
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

  const newVarData = new Float32Array(world.nObs);
  for (let i = 0; i < world.nObs; i += 1) {
    newVarData[i] = varData[world.obsAnnotations[i].__index__];
  }
  return newVarData;
}
