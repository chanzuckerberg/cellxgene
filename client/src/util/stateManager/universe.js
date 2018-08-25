"use strict";
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

    this.schema = {};

    this.obsAnnotations = [];
    this.varAnnotations = [];
    this.obsNameToIndexMap = {}; /* reverse map name to index */
    this.varNameToIndexMap = {}; /* reverse map name to index */

    this.obsLayout = { X: [], Y: [] };
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
  _toSchema(ota) {
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

  _toObsAnnotations(ota) {
    /*
    v0.1 format for metadata:
    metadata: [ { key: val, key: val, ... }, ... ]
    */
    return _.map(ota.data.metadata, (c, i) => {
      return {
        __obsIndex__: i,
        name: c.CellName,
        ...c
      };
    });
  }

  _toVarAnnotations(ota) {
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
    */
    let uz = _.unzip(ota.data.graph);
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
    this.schema = this._toSchema(OTAresponse);
    this.varAnnotations = this._toVarAnnotations(OTAresponse);
    this.init.initialize = true;
    this._tryFinalization();
    return this;
  }

  initFromCells(OTAresponse) {
    /*
    NOTE: this code *assumes* that cell order in data.metadata and data.graph
    are the same.  TODO: error checking.
    */
    this.obsAnnotations = this._toObsAnnotations(OTAresponse);
    this.obsLayout = this._toLayout(OTAresponse);

    this.init.cells = true;
    this._tryFinalization();
    return this;
  }
}

export { UniverseV01 as Universe };
