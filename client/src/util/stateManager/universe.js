// jshint esversion: 6

import _ from "lodash";

/*
This is the public API.   Any non-underscore key in the object is public.
*/
class UniverseBase {
  /*
  creates empty universe, noting which REST API version it presumes
  */
  constructor(apiVersion) {
    this.api = apiVersion;

    /* data status */
    this.finalized = false;
    this.loading = false;
    this.error = null;

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
    Expression dataframes for var/gene. Sparse.  Object keyed by varIndex.
    XXX TODO: periodically remove those not in use so this is a cache.
    */
    this.varData = {
      /* eg, 383: [ 3, 49, 9, ... ] */
    };
  }

  /*
  Universe is built - freeze and generate any derivative information we
  compute on the client side (eg, ranges).
  */
  _finalize() {
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

  varDataByName(name) {
    return this.varData[this.varNameToIndexMap[name]];
  }
}

/*
Universe on top of the REST 0.1 interface
*/
class UniverseV01 extends UniverseBase {
  constructor() {
    super("0.1");

    this.init = {
      initialize: false,
      cells: false
    };
  }

  /*
  Cherry pick from /api/v0.1 response format to make somethign similar
  to the v0.2 schema
  */
  static _toSchema(ota) {
    /*
    keep the 0.1 format for now.  omit CellName to be consistent with
    _toAnnotations()
    */
    return {
      annotations: {
        obs: _(ota.data.schema)
          .omit("CellName")
          .map((val, key) => ({
            name: key,
            ...val
          }))
          .value(),
        var: [{ name: "name", type: "string", variabletype: "categorical" }]
      }
    };
  }

  static _toObsAnnotations(ota) {
    /*
    v0.1 format for metadata:
    metadata: [ { key: val, key: val, ... }, ... ]
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

  static _toLayout(ota) {
    /*
    v0.1 format for the graph is:
    [ [ 'cellname', x, y ], [ 'cellname', x, y, ], ... ]
    */
    const uz = _.unzip(ota.data.graph);
    return {
      X: uz[1],
      Y: uz[2]
    };
  }

  _tryFinalization() {
    if (_.every(this.init)) {
      this._finalize();
    }
  }

  initFromInitialize(OTAresponse) {
    this.schema = UniverseV01._toSchema(OTAresponse);
    this.varAnnotations = UniverseV01._toVarAnnotations(OTAresponse);
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
    this.obsAnnotations = UniverseV01._toObsAnnotations(OTAresponse);
    this.obsLayout = UniverseV01._toLayout(OTAresponse);
    this.nObs = this.obsAnnotations.length;

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
      this.varData[this.varNameToIndexMap[gene]] = data;
    }

    return this;
  }
}

export default UniverseV01;
