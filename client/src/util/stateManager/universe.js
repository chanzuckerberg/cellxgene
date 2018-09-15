// jshint esversion: 6

import _ from "lodash";
import * as kvCache from "./keyvalcache";

/*
This module implements functions that support storage of "Universe",
aka all of the var/obs data and annotations.

These functions are used exclusively by the actions and reducers to
build an internal POJO for use by the rendering components.
*/

/*
Cherry pick from /api/v0.1 response format to make somethign similar
to the v0.2 schema, which we use for internal interfaces.
*/
function OtaRESTv01ToSchema(ota) {
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
      obs: _.map(ota.data.schema, (val, key) => {
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

function OtaRESTv01ToVarAnnotations(ota) {
  /*
  v0.1 initialize response contains 'genes' - names of all genes
  in order.
  */
  return _.map(ota.data.genes, (g, i) => ({ __varIndex__: i, name: g }));
}

function OtaRESTv01ToObsAnnotations(ota) {
  /*
  v0.1 format for metadata:
  metadata: [ { key: val, key: val, ... }, ... ]

  Target format is essentially the same, except the CellName key becomes name.
  */
  return _.map(ota.data.metadata, (c, i) => ({
    __obsIndex__: i,
    name: c.CellName,
    ...c
  }));
}

function OtaRESTv01ToLayout(obsAnnotations, ota) {
  /*
  v0.1 format for the graph is:
  [ [ 'cellname', x, y ], [ 'cellname', x, y, ], ... ]

  NOTE XXX: this code does not assume any particular array ordering in the V0.1
  response.  But for Universe initial load, the layout will be in the same
  order as annotations, so this extra work isn't really necessary.
  */

  const obsAnnotationsByName = _.keyBy(obsAnnotations, "name");
  const { graph } = ota.data;
  const layout = {
    X: new Float32Array(graph.length),
    Y: new Float32Array(graph.length)
  };

  for (let i = 0; i < graph.length; i += 1) {
    const [name, x, y] = graph[i];
    const anno = obsAnnotationsByName[name];
    const idx = anno.__obsIndex__;
    layout.X[idx] = x;
    layout.Y[idx] = y;
  }
  return layout;
}

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

function templateUniverse() {
  /* default universe template */
  const VarDataCacheLowWatermark = 32;
  const VarDataCacheTTLMs = 1000;

  return {
    api: "0.1",
    finalized: true, // XXX: may not be needed

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

export function createUniverseFromRESTv01Response(initOTAResp, cellsOTAResp) {
  /*
  build & return universe from a REST 0.1 /init and /cells response
  */

  const universe = templateUniverse();

  /* extract information from init OTA response */
  universe.schema = OtaRESTv01ToSchema(initOTAResp);
  universe.varAnnotations = OtaRESTv01ToVarAnnotations(initOTAResp);
  universe.nVar = universe.varAnnotations.length;

  /* extract information fron cells OTA response */
  /*
  NOTE: this code *assumes* that cell order in data.metadata and data.graph
  are the same.  TODO: error checking.
  */
  universe.obsAnnotations = OtaRESTv01ToObsAnnotations(cellsOTAResp);
  universe.nObs = universe.obsAnnotations.length;
  universe.obsLayout = OtaRESTv01ToLayout(
    universe.obsAnnotations,
    cellsOTAResp
  );

  return finalize(universe);
}

export function convertExpressionRESTv01ToObject(universe, ota) {
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
  const response = {};
  const { genes, cells } = ota.data;
  for (let idx = 0; idx < genes.length; idx += 1) {
    const gene = genes[idx];
    const data = new Float32Array(universe.nObs);
    for (let c = 0; c < cells.length; c += 1) {
      const obsIndex = universe.obsNameToIndexMap[cells[c].cellname];
      data[obsIndex] = cells[c].e[idx];
    }
    response[gene] = data;
  }
  return response;
}
