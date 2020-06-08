import clip from "../clip";
import { layoutDimensionName, obsAnnoDimensionName } from "../nameCreators";
import * as Dataframe from "../dataframe";
import { isContinuousAnnotation } from "./annotationsHelpers";

/*

World is a subset of universe.  Most code should use world, and should
(generally) not use Universe.  World contains any per-obs or per-var data
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

* unclipped: will contain unclipped variants of all potentiall clipped
  dataframes (obsAnnotations, varData).

*/

function templateWorld() {
  const obsAnnotations = Dataframe.Dataframe.empty();
  const varAnnotations = Dataframe.Dataframe.empty();
  const obsLayout = Dataframe.Dataframe.empty();
  const varData = Dataframe.Dataframe.empty(null, new Dataframe.KeyIndex());
  return {
    /* schema/version related */
    schema: null,
    nObs: 0,
    nVar: 0,
    clipQuantiles: { min: 0, max: 1 },

    /* annotations */
    obsAnnotations,
    varAnnotations,

    /* layout of graph. Dataframe. */
    obsLayout,

    /* Var data columns - subset of all data (may be empty) */
    varData,

    /* unclipped dataframes - subset, but not value clipped */
    unclipped: {
      obsAnnotations,
      varData,
    },
  };
}

function clipDataframe(
  df,
  lowerQuantile,
  upperQuantile,
  quantileF,
  clipPredicate = () => true,
  value = Number.NaN
) {
  /*
  For all columns in the dataframe, clip all values above or below specified 
  quantiles to `value` if clipPredicate returns True for that column (if it
  returns false, skip the column entirely).

  Returns a clipped copy - does not mutate original.

  clipPredicate must have signature: (dataframe, colIndex, colLabel) => boolean
  True signifies that the column should be clipped; false indicates that the
  column should be left intact/unchanged.

  quantileF must have signature:  (label, qval) => number
  */
  if (lowerQuantile < 0) lowerQuantile = 0;
  if (upperQuantile > 1) upperQuantile = 1;
  if (lowerQuantile === 0 && upperQuantile === 1) return df;

  const keys = df.colIndex.labels();
  return df.mapColumns((col, colIdx) => {
    const colLabel = keys[colIdx];
    if (!clipPredicate(df, colIdx, colLabel)) return col;

    const colMin = quantileF(colLabel, lowerQuantile);
    const colMax = quantileF(colLabel, upperQuantile);
    const newCol = clip(col.slice(), colMin, colMax, value);
    return newCol;
  });
}

/*
Create World with contents eq entire universe.   Commonly used to initialize World.
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

  /* save unclipped copies of potentially clipped dataframes */
  world.unclipped = {
    obsAnnotations: world.obsAnnotations.clone(),
    varData: world.varData.clone(),
  };

  return world;
}

/*
clip dataframes based on quantiles.

This is an in-place operation on the world object provided as an argument.
The values in world.unclipped are clipped and assigned to world.obsAnnotations
and world.varData.
*/
function setClippedDataframes(world) {
  const { schema } = world;
  const isContinuousObsAnnotation = (df, idx, label) =>
    isContinuousAnnotation(schema, label);
  const obsQuantile = (label, q) =>
    world.unclipped.obsAnnotations.col(label).summarize().percentiles[100 * q];
  world.obsAnnotations = clipDataframe(
    world.unclipped.obsAnnotations,
    world.clipQuantiles.min,
    world.clipQuantiles.max,
    obsQuantile,
    isContinuousObsAnnotation
  );

  const varDataQuantile = (label, q) =>
    world.unclipped.varData.col(label).summarize().percentiles[100 * q];
  world.varData = clipDataframe(
    world.unclipped.varData,
    world.clipQuantiles.min,
    world.clipQuantiles.max,
    varDataQuantile,
    () => true
  );
}

/*
Subset the current world based upon the current selection, maintaining any existing
clip.  Returns new world.  Parameters:
  * universe
  * world - the current world
  * crossfilter - the selection state
*/
export function createWorldBySelection(universe, world, crossfilter) {
  const newWorld = { ...world, obsLayout: null, unclipped: {}, varData: null };

  /* subset unclipped dataframes based upon current selection */
  const mask = crossfilter.allSelectedMask();
  newWorld.obsLayout = world.obsLayout.isubsetMask(mask);
  newWorld.unclipped.obsAnnotations = world.unclipped.obsAnnotations.isubsetMask(
    mask
  );
  if (world.unclipped.varData.isEmpty()) {
    newWorld.unclipped.varData = world.unclipped.varData.clone();
  } else {
    newWorld.unclipped.varData = world.unclipped.varData.isubsetMask(mask);
  }
  /* subsetting changings dimension size */
  newWorld.nObs = newWorld.unclipped.obsAnnotations.dims[0];

  /* and now clip */
  setClippedDataframes(newWorld);
  return newWorld;
}

/*
Change clip quantiles on the current world, returning a new world.
Parameters:
  * universe
  * world - current world
  * clipQuantiles - new clip
*/
export function createWorldWithNewClip(
  universe,
  world,
  crossfilter,
  clipQuantiles
) {
  const newWorld = { ...world, obsAnnotation: null, varData: null };
  newWorld.clipQuantiles = clipQuantiles;
  newWorld.obsLayout = world.obsLayout.clone();
  newWorld.unclipped = {
    obsAnnotations: world.unclipped.obsAnnotations.clone(),
    varData: world.unclipped.varData.clone(),
  };

  /* and now clip */
  setClippedDataframes(newWorld);
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

function addObsDimension(crossfilter, world, anno) {
  /*
  add single dimension to the crosfilter
  */
  const { obsAnnotations } = world;
  if (obsAnnotations.hasCol(anno.name)) {
    const dimType = deduceDimensionType(anno, anno.name);
    const colData = obsAnnotations.col(anno.name).asArray();
    const name = obsAnnoDimensionName(anno.name);
    if (dimType === "enum") {
      return crossfilter.addDimension(name, "enum", colData);
    }
    if (dimType) {
      return crossfilter.addDimension(name, "scalar", colData, dimType);
    }
  }
  return crossfilter;
}

export function addObsDimensions(crossfilter, world) {
  /*
  Add to crossfilter any dimension present in world.obsAnnotations
  but not yet in the crossfilter
  */
  const schema = world.schema.annotations.obsByName;
  const dimsWeNeed = world.obsAnnotations.colIndex.labels();
  crossfilter = dimsWeNeed.reduce((xfltr, name) => {
    const dimName = obsAnnoDimensionName(name);
    if (xfltr.hasDimension(dimName)) return xfltr;
    return addObsDimension(xfltr, world, schema[name]);
  }, crossfilter);
  return crossfilter;
}

export function createObsDimensions(crossfilter, world, XYdimNames) {
  /*
  create and return a crossfilter with a dimension for every obs annotation
  for which we have a supported type, *except* for the index column, indicated
  by schema.annotations.obs.index.
  */
  const { schema, obsLayout } = world;
  const indexName = schema.annotations.obs.index;
  const annoList = schema.annotations.obs.columns.filter(
    (anno) => anno.name !== indexName
  );
  crossfilter = annoList.reduce((xfltr, anno) => {
    return addObsDimension(xfltr, world, anno);
  }, crossfilter);

  return crossfilter.addDimension(
    layoutDimensionName("XY"),
    "spatial",
    obsLayout.col(XYdimNames[0]).asArray(),
    obsLayout.col(XYdimNames[1]).asArray()
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
  const keys = crossfilter.data.rowIndex.labels(); // row keys, aka universe rowIndex

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
