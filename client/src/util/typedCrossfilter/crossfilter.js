import { polygonContains } from "d3";

import PositiveIntervals from "./positiveIntervals";
import BitArray from "./bitArray";
import {
  sortArray,
  lowerBound,
  lowerBoundIndirect,
  upperBoundIndirect
} from "./sort";
import { makeSortIndex } from "./util";

class NotImplementedError extends Error {
  constructor(...params) {
    super(...params);

    // Maintains proper stack trace for where our error was thrown (only available on V8)
    if (Error.captureStackTrace) {
      Error.captureStackTrace(this, NotImplementedError);
    }
  }
}

export default class ImmutableTypedCrossfilter {
  constructor(data, dimensions = {}, selectionCache = null) {
    /*
    Typically, parameter 'data' is one of:
      - Array of objects/records
      - Dataframe (util/dataframe)
    Other parameters are only used internally.

    Object field description:
    - data: reference to the array of records in the crossfilter
    - selectionBitArray: bit array containing the flatted selection state
      of all dimensions.  This is lazily created and is effectively
      a perfomance cache.  Methods which return a new crossfilter,
      such as select(), addDimention() and delDimension(), will pass
      the cache forward to the new object, as the typical "immutable API"
      usage pattern is to retain the new crossfilter and discard the old.
    - dimensions: contains each dimension and its current state:
        - id: bit offset in the cached bit array
        - dim: the dimension object
        - name: the dimension name
        - selection: the dimension's current selection
    */
    this.data = data;
    this.selectionCache = selectionCache; /* BitArray */
    this.dimensions = dimensions; /* name: { id, dim, name, selection } */
  }

  size() {
    return this.data.length;
  }

  all() {
    return this.data;
  }

  dimensionNames() {
    /* return array of all dimensions (by name) */
    return Object.keys(this.dimensions);
  }

  hasDimension(name) {
    return !!this.dimensions[name];
  }

  addDimension(name, type, ...rest) {
    /*
    Add a new dimension to this crossfilter, of type DimensionType.
    Remainder of parameters are dimension-type-specific.
    */
    const { data, selectionCache } = this;

    if (this.dimensions[name] !== undefined) {
      throw new Error(`Adding duplicate dimension name ${name}`);
    }

    this.selectionCache = null; // pass ownership to new crossfilter

    let id;
    if (selectionCache) {
      id = selectionCache.allocDimension();
      selectionCache.selectAll(id);
    }
    const DimensionType = DimTypes[type];
    const dim = new DimensionType(name, data, ...rest);
    const dimensions = {
      ...this.dimensions,
      [name]: {
        id,
        dim,
        name,
        selection: dim.select({ mode: "all" })
      }
    };
    return new ImmutableTypedCrossfilter(data, dimensions, selectionCache);
  }

  delDimension(name) {
    const { data, selectionCache } = this;
    const dimensions = { ...this.dimensions };
    if (dimensions[name] === undefined) {
      throw new ReferenceError(`Unable to delete unknown dimension ${name}`);
    }

    const { id } = dimensions[name];
    delete dimensions[name];
    this.selectionCache = null; // pass ownership to new crossfilter
    if (selectionCache) {
      selectionCache.freeDimension(id);
    }
    return new ImmutableTypedCrossfilter(data, dimensions, selectionCache);
  }

  select(name, spec) {
    /*
    select on named dimension, as indicated by `spec`.   Spec is an object
    specifying the selection, and must contain at least a `mode` field.

    Examples:
      select("foo", {mode: "all"});
      select("bar", {mode: "none"});
      select("mumble", {mode: "exact", values: "blue"});
      select("mumble", {mode: "exact", values: ["red", "green", "blue"]});
      select("blort", {mode: "range", lo: 0, hi: 999.99});
    */
    const { data, selectionCache } = this;
    this.selectionCache = null;
    const dimensions = { ...this.dimensions };
    const { dim, id, selection: oldSelection } = dimensions[name];
    const newSelection = dim.select(spec);
    newSelection.ranges = PositiveIntervals.canonicalize(newSelection.ranges);
    dimensions[name] = { id, dim, name, selection: newSelection };
    ImmutableTypedCrossfilter._dimSelnHasUpdated(
      selectionCache,
      id,
      newSelection,
      oldSelection
    );
    return new ImmutableTypedCrossfilter(data, dimensions, selectionCache);
  }

  static _dimSelnHasUpdated(selectionCache, id, newSeln, oldSeln) {
    /*
    Selection has updated from oldSeln to newSeln.   Update the
    bit array if it exists.  If not, we will lazy create it when
    needed.
    */
    if (selectionCache) {
      /*
      if both new and old selection use the same index, we can
      perform an incremental update. If the index changed, we have
      to do a suboptimal full deselect/select.
      */
      let adds;
      let dels;
      if (newSeln.index === oldSeln.index) {
        adds = PositiveIntervals.difference(newSeln.ranges, oldSeln.ranges);
        dels = PositiveIntervals.difference(oldSeln.ranges, newSeln.ranges);
      } else {
        // console.log("suboptimal selection update - index changed");
        adds = newSeln.ranges;
        dels = oldSeln.ranges;
      }

      /*
      allow dimensions to return selected ranges in either dimension sort
      order (indirect via index), or in original record order.

      If sort index exists in the dimension, assume sort ordered ranges.
      */
      if (oldSeln.index) {
        dels.forEach(interval =>
          selectionCache.deselectIndirectFromRange(id, oldSeln.index, interval)
        );
      } else {
        dels.forEach(interval =>
          selectionCache.deselectFromRange(id, interval)
        );
      }

      if (newSeln.index) {
        adds.forEach(interval =>
          selectionCache.selectIndirectFromRange(id, newSeln.index, interval)
        );
      } else {
        adds.forEach(interval => selectionCache.selectFromRange(id, interval));
      }
    }
  }

  _getSelectionCache() {
    if (!this.selectionCache) {
      // console.log("...rebuilding crossfilter cache...");
      const selectionCache = new BitArray(this.data.length);
      Object.keys(this.dimensions).forEach(name => {
        const { selection } = this.dimensions[name];
        const id = selectionCache.allocDimension();
        this.dimensions[name].id = id;
        const { ranges, index } = selection;
        ranges.forEach(range => {
          if (index) {
            selectionCache.selectIndirectFromRange(id, index, range);
          } else {
            selectionCache.selectFromRange(id, range);
          }
        });
      });
      this.selectionCache = selectionCache;
    }
    return this.selectionCache;
  }

  allSelected() {
    /*
    return array of all records currently selected by all dimensions
    */
    const selectionCache = this._getSelectionCache();
    const { data } = this;
    if (Array.isArray(data)) {
      const res = [];
      for (let i = 0, len = data.length; i < len; i += 1) {
        if (selectionCache.isSelected(i)) {
          res.push(data[i]);
        }
      }
      return res;
    }
    /* else, Dataframe-like */
    return data.isubsetMask(this.allSelectedMask());
  }

  allSelectedMask() {
    /*
    return Uint8Array containing selection state (truthy/falsey) for each record.
    */
    const selectionCache = this._getSelectionCache();
    return selectionCache.fillBySelection(
      new Uint8Array(this.data.length),
      1,
      0
    );
  }

  countSelected() {
    /*
    return number of records selected on all dimensions
    */
    const selectionCache = this._getSelectionCache();
    return selectionCache.selectionCount();
  }

  isElementSelected(i) {
    /*
    return truthy/falsey if this record is selected on all dimensions
    */
    const selectionCache = this._getSelectionCache();
    return selectionCache.isSelected(i);
  }

  fillByIsSelected(array, selectedValue, deselectedValue) {
    /*
    fill array with one of two values, based upon selection state.
    */
    const selectionCache = this._getSelectionCache();
    return selectionCache.fillBySelection(
      array,
      selectedValue,
      deselectedValue
    );
  }
}

/*
Base dimension object.

A Dimension is an index, accessed via a select() method.  The protocol
for a dimension:
  - constructor - first param is name, remainder is whatever params are
    required to initialize the dimension.
  - select - one and only param is the selection specifier.  Returns an
    array of record IDs.
  - name - the dimension name/label.
*/
class _ImmutableBaseDimension {
  constructor(name) {
    this.name = name;
  }

  select(spec) {
    const { mode } = spec;
    if (mode === undefined) {
      throw new Error("select spec does not contain 'mode'");
    }
    throw new Error(
      `select mode ${mode} not implemented by dimension ${this.name}`
    );
  }
}

class ImmutableScalarDimension extends _ImmutableBaseDimension {
  constructor(name, data, value, ValueArrayType) {
    super(name);

    // Three modes - caller can provide a pre-created value array,
    // a map function which will create it, or another array which
    // will used with an identity map function.
    let array;
    if (value instanceof ValueArrayType) {
      // user has provided the final typed array - just use it
      if (value.length !== data.length) {
        throw new RangeError(
          "ScalarDimension values length must equal crossfilter data record count"
        );
      }
      array = value;
    } else if (value instanceof Function) {
      // Create value array from user-provided map function.
      array = this._createValueArray(
        data,
        value,
        new ValueArrayType(data.length)
      );
    } else if (isArrayOrTypedArray(value)) {
      // Create value array from user-provided array.  Typically used
      // only by enumerated dimensions
      array = this._createValueArray(
        data,
        i => value[i],
        new ValueArrayType(data.length)
      );
    } else {
      throw new NotImplementedError(
        "dimension value must be function or value array type"
      );
    }
    this.value = array;

    // create sort index
    this.index = makeSortIndex(array);
  }

  /* eslint-disable class-methods-use-this */
  _createValueArray(data, mapf, array) {
    // create dimension value array
    const len = data.length;
    const larray = array;
    for (let i = 0; i < len; i += 1) {
      larray[i] = mapf(i, data);
    }
    return larray;
  }
  /* eslint-enable class-methods-use-this */

  select(spec) {
    const { mode } = spec;
    const { index } = this;
    switch (mode) {
      case "all":
        return { ranges: [[0, this.value.length]], index };
      case "none":
        return { ranges: [], index };
      case "exact":
        return this.selectExact(spec);
      case "range":
        return this.selectRange(spec);
      default:
        return super.select(spec);
    }
  }

  selectExact(spec) {
    const { value, index } = this;
    let { values } = spec;
    if (!Array.isArray(values)) {
      values = [values];
    }
    const ranges = [];
    for (let v = 0, len = values.length; v < len; v += 1) {
      const r = [
        lowerBoundIndirect(value, index, values[v], 0, value.length),
        upperBoundIndirect(value, index, values[v], 0, value.length)
      ];
      if (r[0] <= r[1]) {
        ranges.push(r);
      }
    }
    return { ranges, index };
  }

  selectRange(spec) {
    const { value, index } = this;
    /* [lo, hi) */
    const { lo, hi } = spec;
    const ranges = [];
    const r = [
      lowerBoundIndirect(value, index, lo, 0, value.length),
      lowerBoundIndirect(value, index, hi, 0, value.length)
    ];
    if (r[0] < r[1]) ranges.push(r);
    return { ranges, index };
  }
}

class ImmutableEnumDimension extends ImmutableScalarDimension {
  constructor(name, data, value) {
    super(name, data, value, Uint32Array);
  }

  _createValueArray(data, mapf, array) {
    const len = data.length;
    const larray = array;

    // create enumeration table - mapping between the value
    // and the enum.
    const s = new Set();
    for (let i = 0; i < len; i += 1) {
      s.add(mapf(i, data));
    }
    const enumIndex = sortArray(Array.from(s));
    this.enumIndex = enumIndex;

    // create dimension value array
    const enumLen = enumIndex.length;
    for (let i = 0; i < len; i += 1) {
      const v = mapf(i, data);
      const e = lowerBound(enumIndex, v, 0, enumLen);
      larray[i] = e;
    }
    return larray;
  }

  selectExact(spec) {
    const { enumIndex } = this;
    const { values } = spec;
    return super.selectExact({
      mode: spec.mode,
      values: values.map(v => lowerBound(enumIndex, v, 0, enumIndex.length))
    });
  }

  /* eslint-disable class-methods-use-this */
  selectRange() {
    throw new Error("range selection unsupported on Enumerated dimension");
  }
  /* eslint-enable class-methods-use-this */
}

class ImmutableSpatialDimension extends _ImmutableBaseDimension {
  constructor(name, data, X, Y) {
    super(name);

    if (X.length !== Y.length && X.length !== data.length) {
      throw new RangeError(
        "SpatialDimension values must have same dimensionality as crossfilter"
      );
    }
    this.X = X;
    this.Y = Y;

    this.Xindex = makeSortIndex(X);
    this.Yindex = makeSortIndex(Y);
  }

  select(spec) {
    const { mode } = spec;
    switch (mode) {
      case "all":
        return { ranges: [[0, this.X.length]], index: null };
      case "none":
        return { ranges: [], index: null };
      case "within-rect":
        return this.selectWithinRect(spec);
      case "within-polygon":
        return this.selectWithinPolygon(spec);
      default:
        return super.select(spec);
    }
  }

  selectWithinRect(spec) {
    /*
      { mode: "within-rect", minX: 1, minY: 0, maxX: 3, maxY: 9 }
    */
    const { minX, minY, maxX, maxY } = spec;
    const { X, Y } = this;
    const ranges = [];
    let start = -1;
    for (let i = 0, l = X.length; i < l; i += 1) {
      const x = X[i];
      const y = Y[i];
      const inside = minX <= x && x < maxX && minY <= y && y < maxY;
      if (inside && start === -1) start = i;
      if (!inside && start !== -1) {
        ranges.push([start, i]);
        start = -1;
      }
    }
    if (start !== -1) ranges.push([start, X.length]);
    return { ranges, index: null };
  }

  /*
  Relatively brute force filter by polygon.

  Currently uses d3.polygonContains() to test for polygon inclusion, which itself
  uses a ray casting (crossing number) algorithm.  There are a series of optimizations
  to make this faster:
    * first sliced by X or Y, using an index on the axis
    * then the polygon bounding box is used for trivial rejection
    * then the polygon test is applied
  */

  selectWithinPolygon(spec) {
    /*
      { mode: "within-polygon", polygon: [ [x0, y0], ... ] }
    */
    const { polygon } = spec;
    const [minX, minY, maxX, maxY] = polygonBoundingBox(polygon);
    const { X, Y, Xindex, Yindex } = this;
    const { length } = X;
    let slice;
    let index;
    if (maxY - minY > maxX - minX) {
      slice = [
        lowerBoundIndirect(X, Xindex, minX, 0, length),
        lowerBoundIndirect(X, Xindex, maxX, 0, length)
      ];
      index = Xindex;
    } else {
      slice = [
        lowerBoundIndirect(Y, Yindex, minY, 0, length),
        lowerBoundIndirect(Y, Yindex, maxY, 0, length)
      ];
      index = Yindex;
    }

    const ranges = [];
    let start = -1;
    for (let i = slice[0], e = slice[1]; i < e; i += 1) {
      const rid = index[i];
      const x = X[rid];
      const y = Y[rid];
      const inside =
        minX <= x &&
        x < maxX &&
        minY <= y &&
        y < maxY &&
        withinPolygon(polygon, x, y);

      if (inside && start === -1) start = i;
      if (!inside && start !== -1) {
        ranges.push([start, i]);
        start = -1;
      }
    }
    if (start !== -1) ranges.push([start, slice[1]]);
    return { ranges, index };
  }
}

/* Helpers */
export const DimTypes = {
  scalar: ImmutableScalarDimension,
  enum: ImmutableEnumDimension,
  spatial: ImmutableSpatialDimension
};

function isArrayOrTypedArray(x) {
  return (
    Array.isArray(x) ||
    (ArrayBuffer.isView(x) &&
      Object.prototype.toString.call(x) !== "[object DataView]")
  );
}

/* return bounding box of the polygon */
function polygonBoundingBox(polygon) {
  let minX = Number.MAX_VALUE;
  let minY = Number.MAX_VALUE;
  let maxX = Number.MIN_VALUE;
  let maxY = Number.MIN_VALUE;
  for (let i = 0, l = polygon.length; i < l; i += 1) {
    const point = polygon[i];
    const [x, y] = point;
    if (x < minX) minX = x;
    if (y < minY) minY = y;
    if (x > maxX) maxX = x;
    if (y > maxY) maxY = y;
  }
  return [minX, minY, maxX, maxY];
}

function withinPolygon(polygon, x, y) {
  // TODO XXX replace
  return polygonContains(polygon, [x, y]);
}
