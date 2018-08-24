"use strict";
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

import PositiveIntervals from "./positiveIntervals";
import BitArray from "./bitArray";
import {
  fillRange,
  lowerBound,
  lowerBoundIndirect,
  upperBound,
  upperBoundIndirect
} from "./util";

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

  dimension(value, valueArrayType) {
    const id = this.selection.allocDimension();
    let dim;
    if (valueArrayType === "enum") {
      dim = new EnumDimension(value, this, id);
    } else {
      dim = new ScalarDimension(value, valueArrayType, this, id);
    }
    this.filters.push({ id, dim });
    dim.filterAll();
    return dim;
  }

  _freeDimension(id) {
    this.selection.freeDimension(id);
    this.filters = this.filters.filter(f => f._id != id);
  }

  // return array of all records that are selected/filtered
  // by all dimensions.
  allFiltered() {
    const selection = this.selection;
    const res = [];
    for (let i = 0, len = this.data.length; i < len; i++) {
      if (selection.isSelected(i)) {
        res.push(this.data[i]);
      }
    }
    return res;
  }

  countFiltered() {
    return this.selection.selectionCount;
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

// Base dimension type - value must be a scalar type (eg, int, float),
// and value array must be a TypedArray.
//
class ScalarDimension {
  constructor(value, valueArrayType, crossfilter, id) {
    this.crossfilter = crossfilter;
    this._id = id;

    // current selection filter, expressed as PostiveIntervals.
    this.currentFilter = [];

    // Create value array
    const array = this._createValueArray(
      value,
      new valueArrayType(this.crossfilter.data.length)
    );
    this.value = array;

    // create sort index
    this.index = fillRange(new Uint32Array(this.crossfilter.data.length));
    this.index.sort((a, b) => array[a] - array[b]);

    // groups, if any
    this.groups = [];
  }

  _createValueArray(value, array) {
    // create dimension value array
    const data = this.crossfilter.data;
    const len = data.length;
    for (let i = 0; i < len; i++) {
      array[i] = value(data[i]);
    }
    return array;
  }

  dispose() {
    this.crossfilter._freeDimension(this._id);
    return this;
  }

  id() {
    return this._id;
  }

  // Argument is an array of intervals indicating records newly selected/filtered
  //
  _updateFilters(newFilter) {
    newFilter = PositiveIntervals.canonicalize(newFilter);

    const adds = PositiveIntervals.difference(newFilter, this.currentFilter);
    const dels = PositiveIntervals.difference(this.currentFilter, newFilter);

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

    this.currentFilter = newFilter;
    this.crossfilter.updateTime += 1;
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
    for (let v = 0, len = values.length; v < len; v++) {
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
    const data = this.crossfilter.data;
    const selection = this.crossfilter.selection;
    const index = this.index;
    const len = index.length;
    const ret = [];
    let i = 0;
    let skip = 0;
    let found = 0;

    // skip up to offset records
    for (i = len - 1; 0 <= i && skip < offset; i--) {
      if (selection.isSelected(index[i])) {
        skip++;
      }
    }

    // grab up to k records
    for (; 0 <= i && found < k; i--) {
      if (selection.isSelected(index[i])) {
        ret.push(data[index[i]]);
        found++;
      }
    }

    return ret;
  }

  // return bottom k records, starting with offset, in ascending order.
  // Order is this dimension's sort order
  bottom(k, offset = 0) {
    const data = this.crossfilter.data;
    const selection = this.crossfilter.selection;
    const index = this.index;
    const len = index.length;
    const ret = [];
    let skip = 0;
    let found = 0;
    let i = 0;

    // skip up to offset records
    for (i = 0; i < len && skip < offset; i++) {
      if (selection.isSelected(index[i])) {
        skip++;
      }
    }

    // grab up to k records
    for (; i < len && found < k; i++) {
      if (selection.isSelected(index[i])) {
        ret.push(data[index[i]]);
        found++;
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
  constructor(value, crossfilter, id) {
    super(value, Uint32Array, crossfilter, id);
  }

  _createValueArray(value, array) {
    const data = this.crossfilter.data;
    const len = data.length;

    // create enumeration table - mapping between the value
    // and the enum.
    const s = new Set();
    for (let i = 0; i < len; i++) {
      s.add(value(data[i]));
    }
    this.enumIndex = Array.from(s);
    this.enumIndex.sort();

    // create dimension value array
    const enumLen = this.enumIndex.length;
    for (let i = 0; i < len; i++) {
      const v = value(data[i]);
      const e = lowerBound(this.enumIndex, v, 0, enumLen);
      array[i] = e;
    }
    return array;
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

// Groups!  Map/reduce
//
class ScalarGroup {
  constructor(groupValue, groupValueType, dimension) {
    // parent dimension
    this.dimension = dimension;

    // generate group names from dimension values
    this.mapValue = this._map(groupValue, groupValueType, dimension);

    // group index is mapping from data record index to group index
    this.groupIndex = new Uint32Array(dimension.crossfilter.data.length);

    // default to counting
    this.reduceCount();

    // Creates this.groups
    this._reduce();
  }

  // internal support function - map all dimension values to group values.
  //
  _map(groupValue, groupValueType, dimension) {
    // groupValue is optional.  Defaults to identity.  Used to perform
    // initial map operation.
    //
    // identity: save some memory...
    if (groupValue === undefined) return dimension.value;

    const data = dimension.value;
    const len = data.length;
    const mapValue = new groupValueType(dimension.value.length);
    for (let i = 0; i < len; i++) {
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
    const selection = this.dimension.crossfilter.selection;
    const data = this.dimension.crossfilter.data;
    intv.forEach(rng => {
      for (let r = rng[0]; r < rng[1]; r++) {
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
    const selection = this.dimension.crossfilter.selection;
    const data = this.dimension.crossfilter.data;
    intv.forEach(rng => {
      for (let r = rng[0]; r < rng[1]; r++) {
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
    const dimension = this.dimension;
    const data = dimension.crossfilter.data;

    // Create groups
    const groupNames = new Set(this.mapValue);
    this.groups = [];
    const groupIndexByName = {};
    groupNames.forEach(name => {
      this.groups.push({ key: name, value: this.reduceInitial() });
      groupIndexByName[name] = this.groups.length - 1;
    });

    // Create groupIndex - index map between data record index and group index
    for (let i = 0, len = this.mapValue.length; i < len; i++) {
      this.groupIndex[i] = groupIndexByName[this.mapValue[i]];
    }

    // reduce all filtered records, IGNORING the current dimension's filter
    const selection = dimension.crossfilter.selection;
    for (let i = 0, len = data.length; i < len; i++) {
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
    return this.reduce((p, v) => p + 1, (p, v) => p - 1, () => 0);
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
  constructor(groupValue, groupValueType, dimension) {
    super(groupValue, groupValueType, dimension);
  }

  _map(groupValue, groupValueType, dimension) {
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

export default crossfilter;
