// jshint esversion: 6

import _ from "lodash";
import * as kvCache from "./keyvalcache";

/*
Private helper function - create and return a template Universe
*/
function templateUniverse() {
  /* default universe template */
  const VarDataCacheLowWatermark = 32;
  const VarDataCacheTTLMs = 1000;

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

    obsLayout: { X: [], Y: [] } /* xy layout */,

    varDataCache: kvCache.create(
      VarDataCacheLowWatermark,
      VarDataCacheTTLMs
    ) /* cache of var data (expression) */
  };
}

/*
This module implements functions that support storage of "Universe",
aka all of the var/obs data and annotations.

These functions are used exclusively by the actions and reducers to
build an internal POJO for use by the rendering components.
*/

// XXX cleanup
/*
Cherry pick from /api/v0.1 response format to make somethign similar
to the v0.2 schema, which we use for internal interfaces.
*/
function RESTv01ResponseToSchema(response) {
  /*
      Annotation schemas in V02 (our target) look like:

          annotations: {
            obs: [
              { name: "name", type: "string" },
              { name: "num_reads", type: "int32" },
              {
                name: "clusters",
                type: "categorical",
                categories=[ 99, 1, "unknown cluster" ]
              },
              { name: "QScore", type: "float32" }
            ],
            var: [
              { "name": "name", "type": "string" },
              { "name": "gene", "type": "string" }
            ]
          }

    In V01, our source, it looks like:

        "schema": {
          "CellName": {
            "displayname": "Name",
            "include": true,
            "type": "string",
            "variabletype": "categorical"
          },
          "Cluster_2d": {
            "displayname": "Cluster2d",
            "include": true,
            "type": "string",
            "variabletype": "categorical"
          },
          "ERCC_reads": {
            "displayname": "ERCC Reads",
            "include": true,
            "type": "int",
            "variabletype": "continuous"
          },
          ...
        }

    Mapping between the two assumes:
      - V01 only has schema for observations
      - CellName is mapped to 'name'
      - type conversion:   float->float32, int->int32, string->string

    */
  return {
    annotations: {
      obs: _.map(response.data.schema, (val, key) => {
        const name = key === "CellName" ? "name" : key;
        let { type } = val;
        if (type === "int") {
          type = "int32";
        }
        if (type === "float") {
          type = "float32";
        }
        return {
          name,
          type
        };
      }),
      var: [{ name: "name", type: "string" }]
    }
  };
}

// XXX cleanup
function RESTv01ResponseToVarAnnotations(response) {
  /*
  v0.1 initialize response contains 'genes' - names of all genes
  in order.
  */
  return _.map(response.data.genes, (g, i) => ({ __varIndex__: i, name: g }));
}

// XXX cleanup
function RESTv01ResponseToObsAnnotations(response) {
  /*
  v0.1 format for metadata:
  metadata: [ { key: val, key: val, ... }, ... ]

  Target format is essentially the same, except the CellName key becomes name.
  */
  return _.map(response.data.metadata, (c, i) => ({
    __index__: i,
    name: c.CellName,
    ...c
  }));
}

// XXX cleanup
function RESTv01ResponseToLayout(obsAnnotations, response) {
  /*
  v0.1 format for the graph is:
  [ [ 'cellname', x, y ], [ 'cellname', x, y, ], ... ]

  NOTE XXX: this code does not assume any particular array ordering in the V0.1
  response.  But for Universe initial load, the layout will be in the same
  order as annotations, so this extra work isn't really necessary.
  */

  const obsAnnotationsByName = _.keyBy(obsAnnotations, "name");
  const { graph } = response.data;
  const layout = {
    X: new Float32Array(graph.length),
    Y: new Float32Array(graph.length)
  };

  for (let i = 0; i < graph.length; i += 1) {
    const [name, x, y] = graph[i];
    const anno = obsAnnotationsByName[name];
    const idx = anno.__index__;
    layout.X[idx] = x;
    layout.Y[idx] = y;
  }
  return layout;
}

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

// XXX cleanup
export function createUniverseFromRESTv01Response(initResponse, cellsResponse) {
  /*
  build & return universe from a REST 0.1 /init and /cells response
  */

  const universe = templateUniverse();

  /* constants */
  universe.api = "0.1";

  /* extract information from init OTA response */
  universe.schema = RESTv01ResponseToSchema(initResponse);
  universe.varAnnotations = RESTv01ResponseToVarAnnotations(initResponse);
  universe.nVar = universe.varAnnotations.length;

  /* extract information fron cells REST json response */
  /*
  NOTE: this code *assumes* that cell order in data.metadata and data.graph
  are the same.  TODO: error checking.
  */
  universe.obsAnnotations = RESTv01ResponseToObsAnnotations(cellsResponse);
  universe.nObs = universe.obsAnnotations.length;
  universe.obsLayout = RESTv01ResponseToLayout(
    universe.obsAnnotations,
    cellsResponse
  );

  return finalize(universe);
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
    .sortBy("__index__")
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

export function createUniverseFromRestV02Response(
  configResponse,
  schemaResponse,
  annotationsObsResponse,
  // XXX: TODO
  // annotationsVarResponse,
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
  universe.schema = schema.annotations;
  universe.nObs = schema.dataframe.nObs;
  universe.nVar = schema.dataframe.nVar;

  /* annotations */
  universe.obsAnnotations = RESTv02AnnotationsResponseToInternal(
    annotationsObsResponse
  );
  // universe.varAnnotations = RESTv02AnnotationsResponseToInternal(
  //   annotationsVarResponse
  // );

  /* layout */
  // To do
  universe.obsLayout = RESTv02LayoutResponseToInternal(layoutObsResponse);

  return finalize(universe);
}

export function convertExpressionRESTv01ToObject(universe, response) {
  /*
    v0.1 ota looks like:
      {
        genes: [ "name1", "name2", ... ],
        cells: [
          { cellname: 'cell1', e: [ 3, 4, n, x, y, ... ] },
          ...
        ]
      }

    convert expression to a simple Float32Array, and return
    [ [geneName, array], [geneName, array], ... ]
    */
  const result = {};
  const { genes, cells } = response.data;
  for (let idx = 0; idx < genes.length; idx += 1) {
    const gene = genes[idx];
    const data = new Float32Array(universe.nObs);
    for (let c = 0; c < cells.length; c += 1) {
      const obsIndex = universe.obsNameToIndexMap[cells[c].cellname];
      data[obsIndex] = cells[c].e[idx];
    }
    result[gene] = data;
  }
  return result;
}
