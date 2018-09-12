// jshint esversion: 6

import _ from "lodash";
import KeyValCache from "./keyvalcache";

/*
This is the public API.   Any non-underscore key in the object is public.
*/
class UniverseBase {
  static VarDataCacheLowWatermark = 32;

  static VarDataCacheTTLMs = 1000;

  /*
  creates empty universe, noting which REST API version it presumes
  */
  constructor(apiVersion) {
    this.api = apiVersion;

    /* data status */
    this.finalized = false;

    this.nObs = 0;
    this.nVar = 0;
    this.schema = {};

    /*
    Annotations
    */
    this.obsAnnotations = []; /* all obs annotations, by obs index */
    this.varAnnotations = []; /* all var annotations, by var index */
    this.obsNameToIndexMap = {}; /* reverse map 'name' to index */
    this.varNameToIndexMap = {}; /* reverse map 'name' to index */

    this.obsLayout = { X: [], Y: [] };

    /*
    Cache for var/gene data (aka expression data), keyed by varIndex (gene).
    Loaded one var at a time.

    Motivation:  the size of the full data (expression) matrix can be quite
    large (numCells X numGenes), and is currently used in a very constrained
    way:  display of differential expression between two genes/vars.

    This LRU cache contains full-dimension expression data for a set of
    variables/genes - ie, an entire column (1-d array) from the expression
    dataframe.

    These var dimensions will be cached until more than "low watermark"
    dimensions are present, and those dimensions are older than TTL.   See
    futher description in the KeyValCache class.
    */
    this.varDataCache = new KeyValCache(
      UniverseBase.VarDataCacheLowWatermark,
      UniverseBase.VarDataCacheTTLMs
    );
  }

  /* shallow clone Universe - used to properly implement reducers */
  clone() {
    return _.clone(this);
  }

  /*
  Universe is built - freeze and generate any derivative information we
  compute on the client side (eg, ranges).
  */
  _finalize() {
    /* A bit of sanity checking! */
    const { nObs, nVar } = this;
    if (
      nObs !== this.obsAnnotations.length ||
      nObs !== this.obsLayout.X.length ||
      nObs !== this.obsLayout.Y.length ||
      nVar !== this.varAnnotations.length
    ) {
      throw new Error("Universe dimensionality mismatch - failed to load");
    }

    this.obsNameToIndexMap = _.transform(
      this.obsAnnotations,
      (acc, value, idx) => {
        acc[value.name] = idx;
      },
      {}
    );
    this.varNameToIndexMap = _.transform(
      this.varAnnotations,
      (acc, value, idx) => {
        acc[value.name] = idx;
      },
      {}
    );
    this.finalized = true;
  }

  varNameToIndex(name) {
    return this.varNameToIndexMap[name];
  }

  varDataByName(name) {
    return this.varDataCache.get(this.varNameToIndexMap[name]);
  }
}

/*
Universe on top of the REST 0.1 interface
*/
class Universe_REST_API_v01 extends UniverseBase {
  constructor() {
    super("0.1");

    this.init = {
      initialize: false,
      cells: false
    };
  }

  /*
  Cherry pick from /api/v0.1 response format to make somethign similar
  to the v0.2 schema.
  */
  static _toSchema(ota) {
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

  static _toObsAnnotations(ota) {
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

  static _toVarAnnotations(ota) {
    /*
    v0.1 initialize response contains 'genes' - names of all genes
    in order.
    */
    return _.map(ota.data.genes, (g, i) => ({ __varIndex__: i, name: g }));
  }

  _toLayout(ota) {
    /*
    v0.1 format for the graph is:
    [ [ 'cellname', x, y ], [ 'cellname', x, y, ], ... ]

    NOTE XXX: this code does not assume any particular array ordering in the V0.1
    response.  But for Universe initial load, the layout will be in the same
    order as annotations, so this extra work isn't really necessary.
    */

    const { obsAnnotations } = this;
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

  _tryFinalization() {
    if (_.every(this.init)) {
      this._finalize();
    }
  }

  initFromInitialize(OTAresponse) {
    this.schema = Universe_REST_API_v01._toSchema(OTAresponse);
    this.varAnnotations = Universe_REST_API_v01._toVarAnnotations(OTAresponse);
    this.nVar = this.varAnnotations.length;
    this.init.initialize = true;
    this._tryFinalization();
    return this;
  }

  initFromCells(OTAresponse) {
    /*
    NOTE: this code *assumes* that cell order in data.metadata and data.graph
    are the same.  TODO: error checking.
    */
    this.obsAnnotations = Universe_REST_API_v01._toObsAnnotations(OTAresponse);
    this.nObs = this.obsAnnotations.length;
    this.obsLayout = this._toLayout(OTAresponse);

    this.init.cells = true;
    this._tryFinalization();
    return this;
  }

  initFromExpression(ota) {
    /*
    v0.1 ota will look like:
      {
        genes: [ "name1", "name2", ... ],
        cells: [
          { cellname: 'cell1', e: [ 3, 4, n, x, y, ... ] },
          ...
        ]
      }
    */
    const { genes, cells } = ota.data;
    for (let idx = 0; idx < genes.length; idx += 1) {
      const gene = genes[idx];
      const data = new Float32Array(this.nObs);
      for (let c = 0; c < cells.length; c += 1) {
        const obsIndex = this.obsNameToIndexMap[cells[c].cellname];
        data[obsIndex] = cells[c].e[idx];
      }
      this.varDataCache.set(this.varNameToIndexMap[gene], data);
    }

    return this;
  }
}

export default Universe_REST_API_v01;
