// jshint esversion: 6

import _ from "lodash";
import * as kvCache from "./keyvalcache";
import summarizeAnnotations from "./summarizeAnnotations";

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

function RESTv02AnnotationsResponseToInternal(response) {
  /*
  Source per the spec:
  {
    names: [
      'tissue_type', 'sex', 'num_reads', 'clusters'
    ],
    data: [
      [ 0, 'lung', 'F', 39844, 99 ],
      [ 1, 'heart', 'M', 83, 1 ],
      [ 49, 'spleen', null, 2, "unknown cluster" ],
      // [ obsOrVarIndex, value, value, value, value ],
      // ...
    ]
  }

  Internal (target) format:
  [
    { __index__: 0, tissue_type: "lung", sex: "F", ... },
    ...
  ]
  */
  const { names, data } = response;
  const keys = ["__index__", ...names];
  return _(data)
    .map(obs => _.zipObject(keys, obs))
    .value();
}

function RESTv02LayoutResponseToInternal(response) {
  /*
  Source per the spec:
  {
    layout: {
      ndims: 2,
      coordinates: [
        [ 0, 0.284483, 0.983744 ],
        [ 1, 0.038844, 0.739444 ],
        // [ obsOrVarIndex, X_coord, Y_coord ],
        // ...
      ]
    }
  }

  Target (internal) format:
  {
    X: Float32Array(numObs),
    Y: Float32Array(numObs)
  }
  In the same order as obsAnnotations
  */
  const { ndims, coordinates } = response.layout;
  if (ndims !== 2) {
    throw new Error("Unsupported layout dimensionality");
  }

  const layout = {
    X: new Float32Array(coordinates.length),
    Y: new Float32Array(coordinates.length)
  };

  for (let i = 0; i < coordinates.length; i += 1) {
    const [idx, x, y] = coordinates[i];
    layout.X[idx] = x;
    layout.Y[idx] = y;
  }
  return layout;
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
  layoutObsResponse
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
  universe.obsAnnotations = RESTv02AnnotationsResponseToInternal(
    annotationsObsResponse
  );
  universe.varAnnotations = RESTv02AnnotationsResponseToInternal(
    annotationsVarResponse
  );

  /* layout */
  universe.obsLayout = RESTv02LayoutResponseToInternal(layoutObsResponse);

  universe.summary = summarizeAnnotations(
    universe.schema,
    universe.obsAnnotations,
    universe.varAnnotations
  );

  reconcileSchemaCategoriesWithSummary(universe);
  return finalize(universe);
}

export function convertExpressionRESTv02ToObject(universe, response) {
  /*
  /data/obs response looks like:
  {
    var: [ varIndices fetched ],
    obs: [
      [ obsIndex, evalue, ... ],
      ...
    ]
  }

  convert expression toa simple Float32Array, and return
  { geneName: array, geneName: array, ... }
  NOTE: geneName, not varIndex
  */
  const vars = response.var;
  const { obs } = response;
  const result = {};
  // XXX TODO: could this use _.unzip and have less code?
  for (let varIdx = 0; varIdx < vars.length; varIdx += 1) {
    const gene = universe.varAnnotations[vars[varIdx]].name;
    const data = new Float32Array(universe.nObs);
    for (let obsIdx = 0; obsIdx < obs.length; obsIdx += 1) {
      data[obsIdx] = obs[obsIdx][varIdx + 1];
    }
    result[gene] = data;
  }
  return result;
}
