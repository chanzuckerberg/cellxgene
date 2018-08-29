// jshint esversion: 6

import _ from "lodash";
import crossfilter from "../typedCrossfilter";
import { parseRGB } from "../parseRGB";
import KeyValCache from "./keyvalcache";
import * as globals from "../../globals";

/*
World is a subset of universe.   Most code should use world, and should
(generally) not use Universe.   World contains any per-obs or per-var data
that must be consisstent acorss the app when we view/manipulate subsets
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
    * __cellIndex__ renamed to __obsIndex__
    * __x__ and __y__ are now in world.obsLayout
    * __color__ and __colorRBG__ should be moved to controls reducer

* obsLayout:

  obsLayout will return an object containing two arrays, containing X and Y
  coordinates respectively.

  Example: { X: [ 0.33, 0.23, ... ], Y: [ 0.8, 0.777, ... ]}

* obsCrossfilter - a crossfilter object across world.obsAnnotations

* obsDimensionMap - an object mapping annotation names to dimensions on
  the obsCrossfilter

*/

/* WorldBase defines hte public API */
class World {
  /*
  create world which is eq universe.  Universe MUST be fully
  initialized for this to suceed (eg all setters called, finalize()
  called)
  */

  static VarDataCacheLowWatermark = 32;

  static VarDataCacheTTLMs = 1000;

  constructor(universe, color) {
    if (!universe.finalized) {
      throw new Error("World can't be created from an partial Universe");
    }
    let defaultColor = color;
    if (!defaultColor || typeof defaultColor !== "string") {
      defaultColor = globals.defaultCellColor;
    }
    // backpointer to our universe
    this._universe = universe;
    // map from the universe obsIndex to our world offset.
    // undefined/null indicates identity map.
    this._worldObsIndex = null;

    /*
    public interface follows
    */

    /* Schema related */
    this.api = universe.api;
    this.schema = universe.schema;
    this.nObs = universe.nObs;
    this.nVar = universe.nVar;

    /* annotations */
    this.obsAnnotations = universe.obsAnnotations;
    this.varAnnotations = universe.varAnnotations;

    /* layout and display characteristics */
    this.obsLayout = universe.obsLayout;
    this.colorName = new Array(this.nObs).fill(defaultColor);

    /* selection management */
    this.obsCrossfilter = crossfilter(this.obsAnnotations);
    this.obsDimensionMap = this._createObsDimensionMap();

    /* derived data & summaries */
    this.summary = this._summarizeAnnotations();
    this.colorRGB = _.map(this.colorName, c => parseRGB(c));

    /*
    cache of var/expression data.   Keyed by var name.
    */
    this.varDataCache = new KeyValCache(
      World.VarDataCacheLowWatermark,
      World.VarDataCacheTTLMs
    );
  }

  /*
  Factory - create world from this world's currently crossfilter selection.
  If you want to create from Universe, just use the constructor.
  */
  static createFromCrossfilterSelection(world, resetColor = false) {
    const newWorld = world.clone();
    const universe = world._universe;
    const resetColorRGB = resetColor ? parseRGB(resetColor) : null;

    /*
    Subset world from universe based upon world's current selection.  Only those
    fields which are subset by observation selection/filtering need to be updated.
    */
    const numSelected = world.obsCrossfilter.countFiltered();
    /*
    If everything is selected, there is nothing to do.  Just return the
    clone object.
    */
    if (numSelected === world.nObs) {
      return newWorld;
    }

    /*
    Else, create a world which is based upon current selection
    */
    newWorld.nObs = numSelected;
    newWorld.obsAnnotations = new Array(numSelected);
    newWorld.obsLayout = {
      X: new Array(numSelected),
      Y: new Array(numSelected)
    };
    newWorld.colorName = new Array(newWorld.nObs);
    newWorld.colorRGB = new Array(newWorld.nObs);
    newWorld._worldObsIndex = new Array(universe.nObs);

    for (let i = 0, sel = 0; i < world.nObs; i += 1) {
      if (world.obsCrossfilter.isElementFiltered(i)) {
        newWorld.obsAnnotations[sel] = world.obsAnnotations[i];
        newWorld.obsLayout.X[sel] = world.obsLayout.X[i];
        newWorld.obsLayout.Y[sel] = world.obsLayout.Y[i];
        if (!resetColor) {
          newWorld.colorName[sel] = world.colorName[i];
          newWorld.colorRGB[sel] = world.colorRGB[i];
        } else {
          newWorld.colorName[sel] = resetColor;
          newWorld.colorRGB[sel] = resetColorRGB;
        }
        sel += 1;
      }
    }

    // build index to our world offset
    newWorld._worldObsIndex.fill(-1); // default - aka unused
    for (let i = 0; i < newWorld.nObs; i += 1) {
      newWorld._worldObsIndex[newWorld.obsAnnotations[i].__obsIndex__] = i;
    }

    newWorld.obsCrossfilter = crossfilter(newWorld.obsAnnotations);
    newWorld.obsDimensionMap = newWorld._createObsDimensionMap();
    newWorld.summary = newWorld._summarizeAnnotations();

    newWorld.varDataCache = new KeyValCache(
      World.VarDataCacheLowWatermark,
      World.VarDataCacheTTLMs
    );
    return newWorld;
  }

  /* return true of the World is eq Universe */
  _worldEqUniverse() {
    return this._universe.obsAnnotations === this.obsAnnotations;
  }

  get obsSelectionUpdateSeq() {
    return this.obsCrossfilter.updateTime;
  }

  varDataByName(name) {
    let vData = this.varDataCache.get(name);
    if (vData !== undefined) {
      return vData;
    }

    const univVarData = this._universe.varDataByName(name);
    if (this._worldEqUniverse()) {
      vData = univVarData;
    } else {
      vData = new Array(this.nObs);
      for (let i = 0; i < this.nObs; i += 1) {
        vData[i] = univVarData[this.obsAnnotations[i].__obsIndex__];
      }
    }
    this.varDataCache.set(name, vData);
    return vData;
  }

  /* shallow clone World - used to properly implement reducers */
  clone() {
    return _.clone(this);
  }

  /*
  set the colors - optionally can provide parsed RBG values
  */
  setColors(colorName, colorRGB) {
    if (colorName.length !== this.nObs) {
      throw new Error("mismatched length - world not consistent");
    }
    if (colorRGB && colorRGB.length !== this.nObs) {
      throw new Error("mismatched length - world not consistent");
    }

    this.colorName = colorName;
    if (!colorRGB) {
      this.colorRGB = _.map(this.colorName, c => parseRGB(c));
    } else {
      this.colorRGB = colorRGB;
    }
    return this;
  }

  _createObsDimensionMap() {
    /*
    create and return a crossfilter dimension for every obs annotation
    for which we have a supported type.
    */
    const { schema, obsCrossfilter, obsLayout, _worldObsIndex } = this;

    const obsDimensionMap = _.transform(
      schema.annotations.obs,
      (result, anno) => {
        const dimType = World._deduceDimensionType(anno, anno.name);
        if (dimType) {
          result[anno.name] = obsCrossfilter.dimension(
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
    const worldIndex = _worldObsIndex ? idx => _worldObsIndex[idx] : idx => idx;
    obsDimensionMap.x = obsCrossfilter.dimension(
      r => obsLayout.X[worldIndex(r.__obsIndex__)],
      Float32Array
    );
    obsDimensionMap.y = obsCrossfilter.dimension(
      r => obsLayout.Y[worldIndex(r.__obsIndex__)],
      Float32Array
    );

    return obsDimensionMap;
  }

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
  _summarizeAnnotations() {
    /*
    Build and return obs/var summary using any annotation in the schema
    */
    const { obsAnnotations } = this;
    const obsSummary = _(this.schema.annotations.obs)
      .keyBy("name")
      .mapValues(anno => {
        const { name, type } = anno;
        const continuous = type === "int32" || type === "float32";

        if (!continuous) {
          return {
            options: _.countBy(obsAnnotations, name)
          };
        }

        if (continuous) {
          let min = Number.POSITIVE_INFINITY;
          let max = Number.NEGATIVE_INFINITY;
          _.forEach(obsAnnotations, obs => {
            const val = Number(obs[name]);
            min = val < min ? val : min;
            max = val > max ? val : max;
          });
          return { range: { min, max } };
        }

        throw new Error("incomprehensible schema");
      })
      .value();

    const varSummary = {}; // TODO XXX - not currently used, so skip it

    return {
      obs: obsSummary,
      var: varSummary
    };
  }

  /*
   Deduce the correct crossfilter dimension type from a metadata
   schema description.
  */
  static _deduceDimensionType(attributes, fieldName) {
    let dimensionType;
    if (attributes.type === "string") {
      dimensionType = "enum";
    } else if (attributes.type === "int32") {
      dimensionType = Int32Array;
    } else if (attributes.type === "float32") {
      dimensionType = Float32Array;
    } else {
      /*
      Currently not supporting boolean and categorical types.
      */
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

export default World;
