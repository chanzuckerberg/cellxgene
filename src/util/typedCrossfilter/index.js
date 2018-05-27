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

var PositiveIntervals = require("./positiveIntervals");
var BitArray = require("./bitArray");
var Util = require("./util");

class TypedCrossfilter {
  constructor(data) {
    this.data = data;

    // filters: array of { id, dimension }
    this.filters = [];
    this.selection = new BitArray(data.length);
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
    this.filters = this.filters.filter(f => f.id != id);
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
    this.id = id;

    // current selection filter, expressed as PostiveIntervals.
    this.currentFilter = [];

    // Create value array
    const array = this._createValueArray(
      value,
      new valueArrayType(this.crossfilter.data.length)
    );
    this.value = array;

    // create sort index
    this.index = Util.fillRange(new Uint32Array(this.crossfilter.data.length));
    this.index.sort((a, b) => array[a] - array[b]);
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
    this.crossfilter._freeDimension(this.id);
  }

  id() {
    return this.id;
  }

  _updateFilters(newFilter) {
    newFilter = PositiveIntervals.canonicalize(newFilter);

    // special case optimization - select all/none can bypass
    // more complex work and just clobber everything.
    //
    if (newFilter.length === 0) {
      this.crossfilter.selection.deselectAll(this.id);
    } else if (
      newFilter.length === 1 &&
      newFilter[0][0] === 0 &&
      newFilter[0][1] == this.index.length
    ) {
      this.crossfilter.selection.selectAll(this.id);
    } else {
      const adds = PositiveIntervals.difference(newFilter, this.currentFilter);
      const dels = PositiveIntervals.difference(this.currentFilter, newFilter);
      dels.forEach(interval =>
        this.crossfilter.selection.deselectIndirectFromRange(
          this.id,
          this.index,
          interval
        )
      );
      adds.forEach(interval =>
        this.crossfilter.selection.selectIndirectFromRange(
          this.id,
          this.index,
          interval
        )
      );
    }

    this.currentFilter = newFilter;
  }

  // filter by value - exact match
  filterExact(value) {
    const newFilter = [
      Util.lowerBoundIndirect(
        this.value,
        this.index,
        value,
        0,
        this.value.length
      ),
      Util.upperBoundIndirect(
        this.value,
        this.index,
        value,
        0,
        this.value.length
      )
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
        Util.lowerBoundIndirect(
          this.value,
          this.index,
          values[v],
          0,
          this.value.length
        ),
        Util.upperBoundIndirect(
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
      Util.lowerBoundIndirect(
        this.value,
        this.index,
        range[0],
        0,
        this.value.length
      ),
      Util.upperBoundIndirect(
        this.value,
        this.index,
        range[1],
        0,
        this.value.length
      )
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
      const e = Util.lowerBound(this.enumIndex, v, 0, enumLen);
      array[i] = e;
    }
    return array;
  }

  filterExact(value) {
    return super.filterExact(
      Util.lowerBound(this.enumIndex, value, 0, this.enumIndex.length)
    );
  }

  filterEnum(values) {
    return super.filterEnum(
      values.map(v =>
        Util.lowerBound(this.enumIndex, v, 0, this.enumIndex.length)
      )
    );
  }

  filterRange(range) {
    return super.filterEnum(
      range.map(v =>
        Util.lowerBound(this.enumIndex, v, 0, this.enumIndex.length)
      )
    );
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

module.exports = crossfilter;
