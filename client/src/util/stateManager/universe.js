// jshint esversion: 6

import _ from "lodash";

import * as kvCache from "./keyvalcache";
import summarizeAnnotations from "./summarizeAnnotations";
import decodeMatrixFBS from "./matrix";

/*
Private helper function - create and return a template Universe
*/
function templateUniverse() {
  /* default universe template */

  /* varDataCache config - see kvCache for semantics */
  const VarDataCacheLowWatermark = 32; // cache element count
  const VarDataCacheTTLMs = 1000; // min cache time in MS

  return {
    api: null,
    finalized: false, // XXX: may not be needed

    nObs: 0,
    nVar: 0,
    schema: {},

    /*
    Annotations
    */
    obsAnnotations: [] /* all obs annotations, by obs index */,
    varAnnotations: [] /* all var annotations, by var index */,
    obsNameToIndexMap: {} /* reverse map 'name' to index */,
    varNameToIndexMap: {} /* reverse map 'name' to index */,
    summary: null /* derived data summaries XXX: consider exploding in place */,

    obsLayout: { X: [], Y: [] } /* xy layout */,

    /*
    Cache of var data (expression), by var annotation name.   Data can be
    accesses as a POJO, but if you want caching semantics, use the kvCache
    API (eg., kvCache.get(), kvCache.set(), ...), which will maintain the
    LRU semantics.
    */
    varDataCache: kvCache.create(VarDataCacheLowWatermark, VarDataCacheTTLMs)
  };
}

/*
This module implements functions that support storage of "Universe",
aka all of the var/obs data and annotations.

These functions are used exclusively by the actions and reducers to
build an internal POJO for use by the rendering components.
*/

/*
generate any client-side transformations or summarization that
is independent of REST API response formats.
*/
function finalize(universe) {
  /* A bit of sanity checking! */
  const { nObs, nVar } = universe;
  if (
    nObs !== universe.obsAnnotations.length ||
    nObs !== universe.obsLayout.X.length ||
    nObs !== universe.obsLayout.Y.length ||
    nVar !== universe.varAnnotations.length
  ) {
    throw new Error("Universe dimensionality mismatch - failed to load");
  }
  // TODO: add more sanity checks, such as:
  //  - all annotations in the schema
  //  - layout has supported number of dimensions
  //  - ...

  /*
  Create all derived (convenience) data structures.
  */
  universe.obsNameToIndexMap = _.transform(
    universe.obsAnnotations,
    (acc, value, idx) => {
      acc[value.name] = idx;
    },
    {}
  );
  universe.varNameToIndexMap = _.transform(
    universe.varAnnotations,
    (acc, value, idx) => {
      acc[value.name] = idx;
    },
    {}
  );
  universe.finalized = true;
  return universe;
}

function RESTv02AnotationsFBSResponseToInternal(arrayBuffer) {
  /*
  Convert a DataFrame FBS to our internal format -- row-major array of
  observations/cells, stored as an object.  Each obs has a key for each
  annotation, plus __index__ containing its obsIndex.

  Example:
  [
    { __index__: 0, tissue_type: "lung", sex: "F", ... },
    ...
  ]

  XXX TODO: we could make use of the columns in building crossfilter
  dimensions (they have to be recreated).  Future optimization.
  */
  const fbs = decodeMatrixFBS(arrayBuffer);
  const keys = fbs.colIdx;
  const result = Array(fbs.nRows);
  for (let row = 0; row < fbs.nRows; row += 1) {
    const rec = { __index__: row };
    for (let col = 0; col < fbs.nCols; col += 1) {
      rec[keys[col]] = fbs.columns[col][row];
    }
    result[row] = rec;
  }
  return result;
}

function RESTv02LayoutFBSResponseToInternal(arrayBuffer) {
  const fbs = decodeMatrixFBS(arrayBuffer, true);
  return {
    X: fbs.columns[0],
    Y: fbs.columns[1]
  };
}

function reconcileSchemaCategoriesWithSummary(universe) {
  /*
  where we treat types as (essentially) categorical metadata, update
  the schema with data-derived categories (in addition to those in
  the server declared schema).

  For example, boolean defined fields in the schema do not contain
  explicit declaration of categories (nor do string fields).  In these
  cases, add a 'categories' field to the schema so it is accessible.
  */

  _.forEach(universe.schema.annotations.obs, s => {
    if (
      s.type === "string" ||
      s.type === "boolean" ||
      s.type === "categorical"
    ) {
      const categories = _.union(
        _.get(s, "categories", []),
        _.get(universe.summary.obs[s.name], "categories", [])
      );
      s.categories = categories;
    }
  });
}

export function createUniverseFromRestV02Response(
  configResponse,
  schemaResponse,
  annotationsObsResponse,
  annotationsVarResponse,
  layoutFBSResponse
) {
  /*
  build & return universe from a REST 0.2 /config, /schema and /annotations/obs response
  */
  const { schema } = schemaResponse;
  const universe = templateUniverse();

  /* constants */
  universe.api = "0.2";

  /* schema related */
  universe.schema = schema;
  universe.nObs = schema.dataframe.nObs;
  universe.nVar = schema.dataframe.nVar;

  /* annotations */
  universe.obsAnnotations = RESTv02AnotationsFBSResponseToInternal(
    annotationsObsResponse
  );
  universe.varAnnotations = RESTv02AnotationsFBSResponseToInternal(
    annotationsVarResponse
  );

  /* layout */
  universe.obsLayout = RESTv02LayoutFBSResponseToInternal(layoutFBSResponse);

  universe.summary = summarizeAnnotations(
    universe.schema,
    universe.obsAnnotations,
    universe.varAnnotations
  );

  reconcileSchemaCategoriesWithSummary(universe);
  return finalize(universe);
}

export function convertDataFBStoObject(universe, arrayBuffer) {
  /*
  /data/var returns a flatbuffer (FBS) as described by cellxgene/fbs/matrix.fbs

  This routine converts the binary wire encoding into a JS object:

  {
    gene: Float32Array,
    ...
  }
  */
  const fbs = decodeMatrixFBS(arrayBuffer);
  const { colIdx, columns } = fbs;
  const result = {};

  for (let c = 0; c < colIdx.length; c += 1) {
    const gene = universe.varAnnotations[colIdx[c]].name;
    result[gene] = columns[c];
  }
  return result;
}
