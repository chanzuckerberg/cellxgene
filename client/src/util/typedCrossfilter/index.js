// jshint esversion: 6

/*
Typedarray Crossfilter - a re-implementation of a subset of crossfilter, with
time/space optimizations predicated upon the following assumptions:
  - dimensions are uniformly typed, and all values must be of that type
  - dimension values must be a primitive type (int, float, string).  Arrays
    or other complex types not supported.
  - dimension creation requires call-provided type declaration
  - no support for adding/removing data to an existing crossfilter.  If you
    want to do that, you have to create the new crossfilter, using the new
    data, from scratch.

The actual backing store for a dimension is a TypedArray, enabling significant
performance improvements over the original crossfilter.

There are also a handful of new methods, primarily to take advantage of the
performance (eg, crossfilter.fillBySelection)

Helpful documents (this module tries to follow the original API as much
as is feasable):
    https://github.com/square/crossfilter/
    http://square.github.io/crossfilter/

There is also a newer, community supported fork of crossfilter, with a
more complex API.  In a few cases, elements of that API were incorporated.
    https://github.com/square/crossfilter/

*/
// XXX replace
import { polygonContains } from "d3";

import PositiveIntervals from "./positiveIntervals";
import BitArray from "./bitArray";
import {
  makeSortIndex,
  lowerBound,
  lowerBoundIndirect,
  upperBoundIndirect
} from "./util";

function isArrayOrTypedArray(x) {
  return (
    Array.isArray(x) ||
    (ArrayBuffer.isView(x) &&
      Object.prototype.toString.call(x) !== "[object DataView]")
  );
}

class NotImplementedError extends Error {
  constructor(...params) {
    super(...params);

    // Maintains proper stack trace for where our error was thrown (only available on V8)
    if (Error.captureStackTrace) {
      Error.captureStackTrace(this, NotImplementedError);
    }
  }
}

class TypedCrossfilter {
  constructor(data) {
    /*
    Typically, data is one of:
      - Array of objects/records
      - Dataframe (util/dataframe)
    */
    this.data = data;

    // filters: array of { id, dimension }
    this.filters = [];
    this.selection = new BitArray(data.length);
    this.updateTime = 0;
  }

  size() {
    return this.data.length;
  }

  all() {
    return this.data;
  }

  /*
  Create a crossfilter dimension, upon which filtering (subselection) can
  be done.  Each dimension is typed, and has a particular set of filtering
  semantics.
    * ScalarDimension - backed by TypedArray values, supporting filtering
      by value (within a value range, or one or more exact values)
    * EnumDimension - backed by an enumeration (eg, strings, bools), filtering
      by one or more enum categories.
    * SpatialDimension - backed by 2D points, filter by containment within
      various shapes (currently supports within Rectangle and within Polygon).
  Call this method to create a dimension, passing arguments appropriate for
  the dimension constructor.
  */
  dimension(DimensionType, ...rest) {
    const id = this.selection.allocDimension();
    const dim = new DimensionType(this, id, ...rest);
    this.filters.push({ id, dim });
    dim.filterAll();
    return dim;
  }

  _freeDimension(id) {
    this.selection.freeDimension(id);
    this.filters = this.filters.filter(f => f._id !== id);
  }

  // return array of all records that are selected/filtered
  // by all dimensions.
  allFiltered() {
    const { data, selection } = this;
    if (Array.isArray(data)) {
      const res = [];
      for (let i = 0, len = data.length; i < len; i += 1) {
        if (selection.isSelected(i)) {
          res.push(data[i]);
        }
      }
      return res;
    }
    /* else, Dataframe-like */
    return data.isubsetMask(this.allFilteredMask());
  }

  // return Uint8array containing selection state (truthy/falsey) for each record.
  //
  allFilteredMask() {
    return this.selection.fillBySelection(
      new Uint8Array(this.data.length),
      1,
      0
    );
  }

  countFiltered() {
    return this.selection.selectionCount();
  }

  isElementFiltered(i) {
    return this.selection.isSelected(i);
  }

  // fill array with one of two values, based upon selection state
  fillByIsFiltered(array, selectedValue, deselectedValue) {
    return this.selection.fillBySelection(
      array,
      selectedValue,
      deselectedValue
    );
  }
}

// Base dimension type - not exported.
class _Dimension {
  constructor(xfltr, id) {
    this.crossfilter = xfltr;
    this._id = id;
    this.groups = [];
  }

  dispose() {
    this.crossfilter._freeDimension(this._id);
    return this;
  }

  id() {
    return this._id;
  }

  _filterUpdate() {
    this.crossfilter.updateTime += 1;
  }
}

// Scalar dimension type - value must be a scalar type (eg, int, float),
// and value array must be a TypedArray.
//
class ScalarDimension extends _Dimension {
  constructor(xfltr, id, value, ValueArrayType) {
    super(xfltr, id);

    // current selection filter, expressed as PostiveIntervals.
    this.currentFilter = [];

    // Two modes - caller can provide a pre-created value array,
    // or a map function which will create it.
    let array;
    if (value instanceof ValueArrayType) {
      // user has provided the final typed array - just use it
      if (value.length !== this.crossfilter.data.length) {
        throw new RangeError(
          "ScalarDimension values length must equal crossfilter data record count"
        );
      }
      array = value;
    } else if (value instanceof Function) {
      // Create value array from user-provided map function.
      array = this._createValueArray(
        value,
        new ValueArrayType(this.crossfilter.data.length)
      );
    } else if (isArrayOrTypedArray(value)) {
      // Create value array from user-provided array.  Typically used
      // only by enumerated dimensions
      array = this._createValueArray(
        i => value[i],
        new ValueArrayType(this.crossfilter.data.length)
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

  _createValueArray(value, array) {
    // create dimension value array
    const { data } = this.crossfilter;
    const len = data.length;
    const larray = array;
    for (let i = 0; i < len; i += 1) {
      larray[i] = value(i, data);
    }
    return larray;
  }

  // Argument is an array of intervals indicating records newly selected/filtered
  //
  _updateFilters(newFilter) {
    const cNewFilter = PositiveIntervals.canonicalize(newFilter);

    const adds = PositiveIntervals.difference(cNewFilter, this.currentFilter);
    const dels = PositiveIntervals.difference(this.currentFilter, cNewFilter);

    this.crossfilter.filters.forEach(f =>
      f.dim.groups.forEach(grp => grp._updateReduceDel(this, dels))
    );

    dels.forEach(interval =>
      this.crossfilter.selection.deselectIndirectFromRange(
        this._id,
        this.index,
        interval
      )
    );

    adds.forEach(interval =>
      this.crossfilter.selection.selectIndirectFromRange(
        this._id,
        this.index,
        interval
      )
    );

    this.crossfilter.filters.forEach(f =>
      f.dim.groups.forEach(grp => grp._updateReduceAdd(this, adds))
    );

    this.currentFilter = cNewFilter;
    this._filterUpdate();
  }

  // filter by value - exact match
  filterExact(value) {
    const newFilter = [
      lowerBoundIndirect(this.value, this.index, value, 0, this.value.length),
      upperBoundIndirect(this.value, this.index, value, 0, this.value.length)
    ];
    if (newFilter[0] <= newFilter[1]) {
      this._updateFilters([newFilter]);
    } else {
      this._updateFilters([]);
    }
    return this;
  }

  // filter by a set of values, eg. enum.
  filterEnum(values) {
    const newFilter = [];
    for (let v = 0, len = values.length; v < len; v += 1) {
      const intv = [
        lowerBoundIndirect(
          this.value,
          this.index,
          values[v],
          0,
          this.value.length
        ),
        upperBoundIndirect(
          this.value,
          this.index,
          values[v],
          0,
          this.value.length
        )
      ];
      if (intv[0] <= intv[1]) newFilter.push(intv);
    }
    this._updateFilters(newFilter);
    return this;
  }

  // filter by value range [lo, hi)
  // lo: inclusive, hi: exclusive
  filterRange(range) {
    const newFilter = [];
    const intv = [
      lowerBoundIndirect(
        this.value,
        this.index,
        range[0],
        0,
        this.value.length
      ),
      upperBoundIndirect(this.value, this.index, range[1], 0, this.value.length)
    ];
    if (intv[0] < intv[1]) newFilter.push(intv);
    this._updateFilters(newFilter);
    return this;
  }

  // select all - equivalent of selecting all in this dimension
  filterAll() {
    this._updateFilters([[0, this.value.length]]);
    return this;
  }

  // select none
  filterNone() {
    this._updateFilters([]);
  }

  // return top k records, starting with offset, in descending order.
  // Order is this dimension's sort order
  top(k, offset = 0) {
    const { data, selection } = this.crossfilter;
    const { index } = this;
    const len = index.length;
    const ret = [];
    let i = 0;
    let skip = 0;
    let found = 0;

    // skip up to offset records
    for (i = len - 1; i >= 0 && skip < offset; i -= 1) {
      if (selection.isSelected(index[i])) {
        skip += 1;
      }
    }

    // grab up to k records
    for (; i >= 0 && found < k; i -= 1) {
      if (selection.isSelected(index[i])) {
        ret.push(data[index[i]]);
        found += 1;
      }
    }

    return ret;
  }

  // return bottom k records, starting with offset, in ascending order.
  // Order is this dimension's sort order
  bottom(k, offset = 0) {
    const { data, selection } = this.crossfilter;
    const { index } = this;
    const len = index.length;
    const ret = [];
    let skip = 0;
    let found = 0;
    let i = 0;

    // skip up to offset records
    for (i = 0; i < len && skip < offset; i += 1) {
      if (selection.isSelected(index[i])) {
        skip += 1;
      }
    }

    // grab up to k records
    for (; i < len && found < k; i += 1) {
      if (selection.isSelected(index[i])) {
        ret.push(data[index[i]]);
        found += 1;
      }
    }

    return ret;
  }

  group(groupValue) {
    const grp = new ScalarGroup(groupValue, this.value.constructor, this);
    this.groups.push(grp);
    return grp;
  }

  _freeGroup(group) {
    this.groups = this.groups.filter(e => e !== group);
  }
}

// Ordered enumeration - supports any sortable enumerable type, eg,
// strings, which can be mapped into an fixed numeric range [0..n).
//
class EnumDimension extends ScalarDimension {
  constructor(xfltr, id, value) {
    super(xfltr, id, value, Uint32Array);
  }

  _createValueArray(value, array) {
    const { data } = this.crossfilter;
    const len = data.length;
    const larray = array;

    // create enumeration table - mapping between the value
    // and the enum.
    const s = new Set();
    for (let i = 0; i < len; i += 1) {
      s.add(value(i, data));
    }
    this.enumIndex = Array.from(s);
    this.enumIndex.sort();

    // create dimension value array
    const enumLen = this.enumIndex.length;
    for (let i = 0; i < len; i += 1) {
      const v = value(i, data);
      const e = lowerBound(this.enumIndex, v, 0, enumLen);
      larray[i] = e;
    }
    return larray;
  }

  filterExact(value) {
    return super.filterExact(
      lowerBound(this.enumIndex, value, 0, this.enumIndex.length)
    );
  }

  filterEnum(values) {
    return super.filterEnum(
      values.map(v => lowerBound(this.enumIndex, v, 0, this.enumIndex.length))
    );
  }

  filterRange(range) {
    return super.filterEnum(
      range.map(v => lowerBound(this.enumIndex, v, 0, this.enumIndex.length))
    );
  }

  group(groupValue) {
    const grp = new EnumGroup(groupValue, this.value.constructor, this);
    this.groups.push(grp);
    return grp;
  }
}

/*
Super simple 2D spatial dimension, supporting basic "filter within"
operations.
*/
class SpatialDimension extends _Dimension {
  constructor(xfltr, id, X, Y) {
    super(xfltr, id);

    if (X.length !== Y.length && X.length !== this.crossfilter.data.length) {
      throw new RangeError(
        "SpatialDimension values must have same dimensionality as crossfilter"
      );
    }
    this.X = X;
    this.Y = Y;

    this.Xindex = makeSortIndex(X);
    this.Yindex = makeSortIndex(Y);
  }

  filterAll() {
    this.crossfilter.selection.selectAll(this._id);
    this._filterUpdate();
  }

  filterNone() {
    this.crossfilter.selection.deselectAll(this._id);
    this._filterUpdate();
  }

  /*
  this could be smarter, but we don't currently use it...
  */
  filterWithinRect(northwest, southeast) {
    const [x0, y0] = northwest;
    const [x1, y1] = southeast;
    const { X, Y } = this;
    const seln = this.crossfilter.selection;
    const { _id } = this;
    seln.deselectAll(_id);
    for (let i = 0, l = this.X.length; i < l; i += 1) {
      const x = X[i];
      const y = Y[i];
      if (x0 <= x && x < x1 && y0 <= y && y < y1) {
        seln.selectOne(_id, i);
      }
    }
    this._filterUpdate();
  }

  /*
  Relatively brute force filter by polygon.   Polygon is array of points, where
  each point is [x,y].   Eg, [[x0,y0], [x1,y1], ...].

  Currently uses d3.polygonContains() to test for polygon inclusion, which itself
  uses a ray casting (crossing number) algorithm.  There are a series of optimizations
  to make this faster:
    * first sliced by X or Y, using an index on the axis
    * then the polygon bounding box is used for trivial rejection
    * then the polygon test is applied
  */
  filterWithinPolygon(polygon) {
    /* return bounding box of the polygon */
    function polygonBoundingBox(pg) {
      let minX = Number.MAX_VALUE;
      let minY = Number.MAX_VALUE;
      let maxX = Number.MIN_VALUE;
      let maxY = Number.MIN_VALUE;
      for (let i = 0, l = pg.length; i < l; i += 1) {
        const p = pg[i];
        const x = p[0];
        const y = p[1];
        if (x < minX) minX = x;
        if (y < minY) minY = y;
        if (x > maxX) maxX = x;
        if (y > maxY) maxY = y;
      }
      return [minX, minY, maxX, maxY];
    }

    const [minX, minY, maxX, maxY] = polygonBoundingBox(polygon);
    const { X, Y } = this;
    let slice;
    let index;
    if (maxY - minY > maxX - minX) {
      slice = [
        lowerBoundIndirect(X, this.Xindex, minX, 0, X.length),
        upperBoundIndirect(X, this.Xindex, maxX, 0, X.length)
      ];
      index = this.Xindex;
    } else {
      slice = [
        lowerBoundIndirect(Y, this.Yindex, minY, 0, Y.length),
        upperBoundIndirect(Y, this.Yindex, maxY, 0, Y.length)
      ];
      index = this.Yindex;
    }

    const seln = this.crossfilter.selection;
    const { _id } = this;
    const testWithin = polygonContains; // d3.polygonContains()
    seln.deselectAll(_id);

    for (let i = slice[0], e = slice[1]; i < e; i += 1) {
      const rid = index[i];
      const x = X[rid];
      const y = Y[rid];
      if (
        minX <= x &&
        x < maxX &&
        minY <= y &&
        y < maxY &&
        testWithin(polygon, [x, y])
      ) {
        seln.selectOne(_id, rid);
      }
    }
    this._filterUpdate();
  }
}

// Groups!  Map/reduce
//
class ScalarGroup {
  constructor(groupValue, groupValueType, dimension) {
    // parent dimension
    this.dimension = dimension;

    // generate group names from dimension values
    this.mapValue = this.constructor._map(
      groupValue,
      groupValueType,
      dimension
    );

    // group index is mapping from data record index to group index
    this.groupIndex = new Uint32Array(dimension.crossfilter.data.length);

    // default to counting
    this.reduceCount();

    // Creates this.groups
    this._reduce();
  }

  // internal support function - map all dimension values to group values.
  //
  static _map(groupValue, GroupValueType, dimension) {
    // groupValue is optional.  Defaults to identity.  Used to perform
    // initial map operation.
    //
    // identity: save some memory...
    if (groupValue === undefined) return dimension.value;

    const data = dimension.value;
    const len = data.length;
    const mapValue = new GroupValueType(dimension.value.length);
    for (let i = 0; i < len; i += 1) {
      mapValue[i] = groupValue(data[i]);
    }
    return mapValue;
  }

  // Update the group reduction incrementally.  Called when *any* dimension filter
  // changes.  Guaranteed to be called AFTER the crossfilter is updated.
  //
  // Arguments:
  //  * dim: the dimension that is changing
  //  * intv: interval list of newly selected values on `dim` (adds)
  //
  _updateReduceAdd(dim, intv) {
    // ignore updates to self, as we don't reduce inclusive of our filter
    if (dim === this.dimension || intv.length === 0) return;

    // Each item in the range was just added to `dim`.  It was NOT previously
    // selected - reduceAdd if it is now selected.
    const { data, selection } = this.dimension.crossfilter;
    intv.forEach(rng => {
      for (let r = rng[0]; r < rng[1]; r += 1) {
        const i = dim.index[r];
        if (selection.isSelectedIgnoringDim(i, this.dimension.id())) {
          const group = this.groups[this.groupIndex[i]];
          group.value = this.reduceAdd(group.value, data[i]);
        }
      }
    });
  }

  // Update the group reduction incrementally.  Called when *any* dimension filter
  // changes.   Guaranteed to be called BEFORE the crossfilter is updated.
  //
  // Arguments:
  //  * dim: the dimension that is changing
  //  * intv: interval list of previously selected values on `dim` (dels)
  //
  _updateReduceDel(dim, intv) {
    // ignore updates to self, as we don't reduce inclusive of our filter
    if (dim === this.dimension || intv.length === 0) return;

    // Each item in the range will be remved from `dim`.  reduceRemove if it
    // is currently selected.
    const { data, selection } = this.dimension.crossfilter;
    intv.forEach(rng => {
      for (let r = rng[0]; r < rng[1]; r += 1) {
        const i = dim.index[r];
        if (selection.isSelectedIgnoringDim(i, this.dimension.id())) {
          const group = this.groups[this.groupIndex[i]];
          group.value = this.reduceRemove(group.value, data[i]);
        }
      }
    });
  }

  // Reduce the entire data set, creating both the group index and the
  // groups data.
  //
  _reduce() {
    const { dimension } = this;
    const { data } = dimension.crossfilter;

    // Create groups
    const groupNames = new Set(this.mapValue);
    this.groups = [];
    const groupIndexByName = {};
    groupNames.forEach(name => {
      this.groups.push({ key: name, value: this.reduceInitial() });
      groupIndexByName[name] = this.groups.length - 1;
    });

    // Create groupIndex - index map between data record index and group index
    for (let i = 0, len = this.mapValue.length; i < len; i += 1) {
      this.groupIndex[i] = groupIndexByName[this.mapValue[i]];
    }

    // reduce all filtered records, IGNORING the current dimension's filter
    const { selection } = dimension.crossfilter;
    for (let i = 0, len = data.length; i < len; i += 1) {
      if (selection.isSelectedIgnoringDim(i, dimension.id())) {
        const group = this.groups[this.groupIndex[i]];
        group.value = this.reduceAdd(group.value, data[i]);
      }
    }
  }

  dispose() {
    this.dimension._freeGroup(this);
    return this;
  }

  // return number of distinct values in the group, independent of any filters.
  //
  size() {
    return this.groups.length;
  }

  // Set the reduce functions and return the grouping.
  //
  reduce(add, remove, initial) {
    this.reduceAdd = add;
    this.reduceRemove = remove;
    this.reduceInitial = initial;
    this._reduce();
    return this;
  }

  // set the reduce functions to count records.
  reduceCount() {
    return this.reduce(p => p + 1, p => p - 1, () => 0);
  }

  // set the reduce functions to sum records using specified value accessor.
  //
  reduceSum(value) {
    return this.reduce((p, v) => p + value(v), (p, v) => p - value(v), () => 0);
  }

  all() {
    const res = [...this.groups];
    res.sort((a, b) => (a.key < b.key ? -1 : a.key > b.key ? 1 : 0));
    return res;
  }
}

class EnumGroup extends ScalarGroup {
  static _map(groupValue, groupValueType, dimension) {
    // groupValue is optional.  Defaults to identity.  Used to perform
    // initial map operation.
    //
    // identity: save some memory
    if (groupValue === undefined) return dimension.value;

    // non-identity mapping unsupported for EnumDimension/EnumGroup.
    // XXX: this could be implemented, but would require another index
    // array to map from the group names/keys back to the dimension values.
    // With this, we just rely on the dimensions `enumIndex` to map from
    // enumeration value to the record.
    throw new NotImplementedError("enumerated group mapping not implemented");
  }

  all() {
    const res = [];
    this.groups.forEach(e =>
      res.push({
        // XXX: assumes identity group map - see comment in _map()
        key: this.dimension.enumIndex[e.key],
        value: e.value
      })
    );
    res.sort((a, b) => (a.key < b.key ? -1 : a.key > b.key ? 1 : 0));
    return res;
  }
}

// Wrapper for backwards compat with crossfilter.
//
function crossfilter(data) {
  return new TypedCrossfilter(data);
}

crossfilter.PositiveIntervals = PositiveIntervals;
crossfilter.BitArray = BitArray;
crossfilter.TypedCrossfilter = TypedCrossfilter;
crossfilter.ScalarDimension = ScalarDimension;
crossfilter.EnumDimension = EnumDimension;
crossfilter.SpatialDimension = SpatialDimension;

export default crossfilter;
