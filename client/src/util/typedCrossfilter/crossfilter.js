// eslint-disable-next-line max-classes-per-file -- classes are interrelated
import PositiveIntervals from "./positiveIntervals";
import BitArray from "./bitArray";
import {
  sortArray,
  lowerBound,
  binarySearch,
  lowerBoundIndirect,
  upperBoundIndirect,
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
  constructor(data, dimensions = {}, selectionCache = {}) {
    /*
    Typically, parameter 'data' is one of:
      - Array of objects/records
      - Dataframe (util/dataframe)
    Other parameters are only used internally.

    Object field description:
    - data: reference to the array of records in the crossfilter
    - selectionCache: object which may contains a bit array indicating
      the flatted selection state for all dimensions, plus other state
      summarizing the selection.  This is lazily created and is effectively
      a perfomance cache.  Methods which return a new crossfilter,
      such as select(), addDimention() and delDimension(), will pass
      the cache forward to the new object, as the typical "immutable API"
      usage pattern is to retain the new crossfilter and discard the old.
      If the cache object is empty, it will be rebuilt.
    - dimensions: contains each dimension and its current state:
        - id: bit offset in the cached bit array
        - dim: the dimension object
        - name: the dimension name
        - selection: the dimension's current selection
    */
    this.data = data;
    this.selectionCache = selectionCache; /* { BitArray, ... }*/
    this.dimensions = dimensions; /* name: { id, dim, name, selection } */
    Object.preventExtensions(this);
  }

  size() {
    return this.data.length;
  }

  all() {
    return this.data;
  }

  setData(data) {
    if (this.data === data) return this;
    // please leave, WIP
    // console.log("...crossfilter set data, will drop cache");
    return new ImmutableTypedCrossfilter(data, this.dimensions);
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
    const { data } = this;
    const { bitArray } = this.selectionCache;

    if (this.dimensions[name] !== undefined) {
      throw new Error(`Adding duplicate dimension name ${name}`);
    }

    this._clearSelectionCache();

    let id;
    if (bitArray) {
      id = bitArray.allocDimension();
      bitArray.selectAll(id);
    }
    const DimensionType = DimTypes[type];
    const dim = new DimensionType(name, data, ...rest);
    Object.freeze(dim);
    const dimensions = {
      ...this.dimensions,
      [name]: {
        id,
        dim,
        name,
        selection: dim.select({ mode: "all" }),
      },
    };

    return new ImmutableTypedCrossfilter(data, dimensions, {
      bitArray,
    });
  }

  delDimension(name) {
    const { data } = this;
    const { bitArray } = this.selectionCache;
    const dimensions = { ...this.dimensions };
    if (dimensions[name] === undefined) {
      throw new ReferenceError(`Unable to delete unknown dimension ${name}`);
    }

    const { id } = dimensions[name];
    delete dimensions[name];
    this._clearSelectionCache();
    if (bitArray) {
      bitArray.freeDimension(id);
    }

    return new ImmutableTypedCrossfilter(data, dimensions, {
      bitArray,
    });
  }

  renameDimension(oldName, newName) {
    const { [oldName]: dim, ...dimensions } = this.dimensions;
    const { data, selectionCache } = this;
    const newDimensions = {
      ...dimensions,
      [newName]: {
        ...dim,
        name: newName,
        dim: dim.dim.rename(newName),
      },
    };
    return new ImmutableTypedCrossfilter(data, newDimensions, selectionCache);
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
    this.selectionCache = {};
    const dimensions = { ...this.dimensions };
    const { dim, id, selection: oldSelection } = dimensions[name];
    const newSelection = dim.select(spec);
    newSelection.ranges = PositiveIntervals.canonicalize(newSelection.ranges);
    dimensions[name] = { id, dim, name, selection: newSelection };
    const newSelectionCache = ImmutableTypedCrossfilter._dimSelnHasUpdated(
      selectionCache,
      id,
      newSelection,
      oldSelection
    );
    return new ImmutableTypedCrossfilter(data, dimensions, newSelectionCache);
  }

  static _dimSelnHasUpdated(selectionCache, id, newSeln, oldSeln) {
    /*
    Selection has updated from oldSeln to newSeln.   Update the
    bit array if it exists.  If not, we will lazy create it when
    needed.
    */
    if (!selectionCache || !selectionCache.bitArray) return {};

    const { bitArray } = selectionCache;

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
      // please leave, WIP
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
      dels.forEach((interval) =>
        bitArray.deselectIndirectFromRange(id, oldSeln.index, interval)
      );
    } else {
      dels.forEach((interval) => bitArray.deselectFromRange(id, interval));
    }

    if (newSeln.index) {
      adds.forEach((interval) =>
        bitArray.selectIndirectFromRange(id, newSeln.index, interval)
      );
    } else {
      adds.forEach((interval) => bitArray.selectFromRange(id, interval));
    }

    return { bitArray };
  }

  _getSelectionCache() {
    if (!this.selectionCache) this.selectionCache = {};

    if (!this.selectionCache.bitArray) {
      // console.log("...rebuilding crossfilter cache...");
      const bitArray = new BitArray(this.data.length);
      Object.keys(this.dimensions).forEach((name) => {
        const { selection } = this.dimensions[name];
        const id = bitArray.allocDimension();
        this.dimensions[name].id = id;
        const { ranges, index } = selection;
        ranges.forEach((range) => {
          if (index) {
            bitArray.selectIndirectFromRange(id, index, range);
          } else {
            bitArray.selectFromRange(id, range);
          }
        });
      });
      this.selectionCache.bitArray = bitArray;
    }
    return this.selectionCache;
  }

  _clearSelectionCache() {
    this.selectionCache = {};
    return this.selectionCache;
  }

  _setSelectionCache(vals = {}) {
    Object.assign(this.selectionCache, vals);
    return this.selectionCache;
  }

  allSelected() {
    /*
    return array of all records currently selected by all dimensions
    */
    const selectionCache = this._getSelectionCache();
    const { bitArray } = selectionCache;
    const { data } = this;
    if (Array.isArray(data)) {
      const res = [];
      for (let i = 0, len = data.length; i < len; i += 1) {
        if (bitArray.isSelected(i)) {
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
    let { allSelectedMask } = selectionCache;

    if (allSelectedMask !== undefined) return allSelectedMask;

    allSelectedMask = selectionCache.bitArray.fillBySelection(
      new Uint8Array(this.data.length),
      1,
      0
    );
    this._setSelectionCache({ allSelectedMask });
    return allSelectedMask;
  }

  countSelected() {
    /*
    return number of records selected on all dimensions
    */
    const selectionCache = this._getSelectionCache();
    let { countSelected } = selectionCache;

    if (countSelected !== undefined) return countSelected;

    countSelected = selectionCache.bitArray.selectionCount();
    this._setSelectionCache({ countSelected });
    return countSelected;
  }

  isElementSelected(i) {
    /*
    return truthy/falsey if this record is selected on all dimensions
    */
    const selectionCache = this._getSelectionCache();
    return selectionCache.bitArray.isSelected(i);
  }

  fillByIsSelected(array, selectedValue, deselectedValue) {
    /*
    fill array with one of two values, based upon selection state.
    */
    const selectionCache = this._getSelectionCache();
    return selectionCache.bitArray.fillBySelection(
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

  clone() {
    return Object.assign(Object.create(Object.getPrototypeOf(this)), this);
  }

  rename(name) {
    const d = this.clone();
    d.name = name;
    return d;
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
        (i) => value[i],
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

  // eslint-disable-next-line class-methods-use-this -- needed for polymorphism
  _createValueArray(data, mapf, array) {
    // create dimension value array
    const len = data.length;
    const larray = array;
    for (let i = 0; i < len; i += 1) {
      larray[i] = mapf(i, data);
    }
    return larray;
  }

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
        upperBoundIndirect(value, index, values[v], 0, value.length),
      ];
      if (r[0] <= r[1]) {
        ranges.push(r);
      }
    }
    return { ranges, index };
  }

  selectRange(spec) {
    const { value, index } = this;
    /* 
    if !inclusive: [lo, hi) else [lo, hi]
    */
    const { lo, hi, inclusive } = spec;
    const ranges = [];
    const r = [
      lowerBoundIndirect(value, index, lo, 0, value.length),
      inclusive
        ? upperBoundIndirect(value, index, hi, 0, value.length)
        : lowerBoundIndirect(value, index, hi, 0, value.length),
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
    let { values } = spec;
    if (!Array.isArray(values)) {
      values = [values];
    }
    return super.selectExact({
      mode: spec.mode,
      values: values.map((v) =>
        binarySearch(enumIndex, v, 0, enumIndex.length)
      ),
    });
  }

  // eslint-disable-next-line class-methods-use-this -- enables polymorphism
  selectRange() {
    throw new Error("range selection unsupported on Enumerated dimension");
  }
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
        lowerBoundIndirect(X, Xindex, maxX, 0, length),
      ];
      index = Xindex;
    } else {
      slice = [
        lowerBoundIndirect(Y, Yindex, minY, 0, length),
        lowerBoundIndirect(Y, Yindex, maxY, 0, length),
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
  spatial: ImmutableSpatialDimension,
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

/**
 *  withinPolygon determines if a point is within a polygon
 *  Code adapted from https://github.com/d3/d3-polygon/blob/master/src/contains.js
 *  @param {array} polygon - is an array of point arrays of format [[x1, y1], [x2, y2], ...]
 *  @param {float} x - point x coordinate
 *  @param {float} y - point y coordinate
 *  @type {boolean}
 */
function withinPolygon(polygon, x, y) {
  const n = polygon.length;
  let p = polygon[n - 1];
  let x0 = p[0];
  let y0 = p[1];
  let x1;
  let y1;
  let inside = false;

  for (let i = 0; i < n; i += 1) {
    p = polygon[i];
    x1 = p[0];
    y1 = p[1];

    if (y1 > y !== y0 > y && x < ((x0 - x1) * (y - y1)) / (y0 - y1) + x1)
      inside = !inside;
    x0 = x1;
    y0 = y1;
  }
  return inside;
}
