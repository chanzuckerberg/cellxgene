"use strict";
// jshint esversion: 6

import _ from "lodash";
import crossfilter from "../typedCrossfilter";

/*
World is a subset of universe.   Most code should use world, and should
(generally) not use Universe.

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
    * __cellIndex__ renamed to __cellId__  XXX: should be in mapToUniverse
    * __x__ and __y__ are now in world.obsLayout
    * __color__ and __colorRBG__ should be moved to controls reducer

* obsLayout:

  obsLayout will return an object containing two arrays, containing X and Y
  coordinates respectively.

  Example: { X: [ 0.33, 0.23, ... ], Y: [ 0.8, 0.777, ... ]}

* matrix:  TBD

* obsCrossfilter - a crossfilter object across world.obsAnnotations

* obsDimensionMap - an object mapping annotation names to dimensions on
  the obsCrossfilter

*/

/* WorldBase defines hte public API */
class WorldBase {
  /*
  create world which is eq universe.  Universe MUST be fully
  initialized for this to suceed (eg all setters called, finalize()
  called)
  */
  constructor(universe) {
    if (!universe.finalized)
      throw new Error("World can't be created from an partial Universe");
    this.api = universe.api;

    this._universe = universe;
    /*
    map from world index to universe index, for cases where world < universe.
    null value signifies world == universe.
    */
    this._universeObsIndex = null;

    /*
    public interface
    */
    this.schema = this._universe.schema;
    this.obsAnnotations = universe.obsAnnotations;
    this.varAnnotations = universe.varAnnotations;
    this.obsLayout = universe.obsLayout;
    this.matrix = undefined; // TODO
    this.obsCrossfilter = crossfilter(this.obsAnnotations);
    this.obsDimensionMap = {};
    this.summary = { obs: [], var: [] };
  }

  get obsSelectionUpdateSeq() {
    return this.obsCrossfilter.updateTime;
  }

  varDataByName(name) {
    return this._universe.varDataByName(name);
  }

  /* Factory - create world from this world's universe */
  createFromUniverse() {
    throw new Error("unimplemented - subclass responsibility");
  }

  /* Factory - create world from this world's currently crossfilter selection */
  createFromSelected() {
    throw new Error("unimplemented - subclass responsibility");
  }
}

class WorldV01 extends WorldBase {
  constructor(universe) {
    super(universe);

    if (this.api != "0.1") throw new Error("unsupported REST API version");
    this._createDimensionMapV01();
    this._summarizeWorldV01();
  }

  _createDimensionMapV01() {
    /*
    create a crossfilter dimension for every obs annotation
    for which we have a supported type.
    */
    this.obsDimensionMap = _.transform(
      this.schema.annotations.obs,
      (result, anno) => {
        const dimType = this._deduceDimensionTypeV01(anno, anno.name);
        if (dimType) {
          result[anno.name] = this.obsCrossfilter.dimension(
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
    this.obsDimensionMap.x = this.obsCrossfilter.dimension(
      r => this.obsLayout.X[r.__obsIndex__],
      Float32Array
    );
    this.obsDimensionMap.y = this.obsCrossfilter.dimension(
      r => this.obsLayout.Y[r.__obsIndex__],
      Float32Array
    );
  }

  /* Factory - create world from this world's universe */
  createFromUniverse() {
    return new WorldV01(this._universe);
  }

  /* Factory - create world from this world's currently crossfilter selection */
  createFromSelected() {
    let world = new WorldV01(this._universe);
    throw new Error("unimplemented");
    // return world;
  }

  /**
   ** Getters for the state inside world.    These are guarnateed to return
   ** data for just the subset included in world.
   **/

  /*
  Summary information for each annotation, keyed by annotation name.
  Value will be an object, containing either 'range' or 'options' object,
  depending on the annotation schema type (categorical or continuous).

  Summarize for BOTH obs and var annotations.  Result format:

  {
    obs: {
      annotation_name: { ... },
      ...
    },
    var: {
      annotation_name: { ... },
      ...
    }
  }

  Example:
    {
      "Splice_sites_Annotated": {
        "range": {
          "min": 26,
          "max": 1075869
        }
      },
      "Selection": {
        "options": {
          "Astrocytes(HEPACAM)": 714,
          "Endothelial(BSC)": 123,
          "Oligodendrocytes(GC)": 294,
          "Neurons(Thy1)": 685,
          "Microglia(CD45)": 1108,
          "Unpanned": 665
        }
      }
    }
  */
  _summarizeWorldV01() {
    // TODO: reminder to future self
    if (this.api != "0.1") throw new Error("unimplemented API version");

    /*
    Build obs summary using any annotation in the obs schema
    */
    const obsAnnotations = this.obsAnnotations;
    const obsSummary = _(this.schema.annotations.obs)
      .keyBy("name")
      .mapValues(anno => {
        const key = anno.name;
        const type = anno.type;
        const continuous = anno.variabletype === "continuous";

        if (!continuous) {
          return {
            options: _.countBy(obsAnnotations, key)
          };
        } else if (continuous && (type === "int" || type === "float")) {
          let min = Number.POSITIVE_INFINITY;
          let max = Number.NEGATIVE_INFINITY;
          _.forEach(obsAnnotations, obs => {
            const val = Number(obs[key]);
            min = val < min ? val : min;
            max = val > max ? val : max;
          });
          return { range: { min, max } };
        } else {
          throw new Error("incomprehensible schema");
        }
      })
      .value();

    /*
    In the V1 REST API, we only have "gene names", a la the 'name' annotation
    */
    const varSummary = {}; // TODO XXX

    this.summary = {
      obs: obsSummary,
      ["var"]: varSummary
    };
  }

  /*
   Deduce the correct crossfilter dimension type from a metadata
   schema description.
  */
  _deduceDimensionTypeV01(attributes, fieldName) {
    let dimensionType;
    if (attributes.type === "string") {
      dimensionType = "enum";
    } else if (attributes.type === "int") {
      dimensionType = Int32Array;
    } else if (attributes.type === "float") {
      dimensionType = Float32Array;
    } else {
      console.error(
        `Warning - REST API returned unknown metadata schema (${
          attributes.type
        }) for field ${fieldName}.`
      );
      // skip it - we don't know what to do with this type
    }
    return dimensionType;
  }
}

export { WorldV01 as World };
