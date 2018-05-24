"use strict";
// jshint esversion: 6

/*
Typedarray Crossfilter - a re-implementation of a subset of crossfilter, with
major time/space optimizations predicated upon the following assumptions:
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

Helpful documents (this code tries to follow the original API as much
as is feasable):
    https://github.com/square/crossfilter/
    http://square.github.io/crossfilter/

There is also a newer, community supported fork of crossfilter, with a
more complex API.  In a few cases, elements of that API were incorporated.
    https://github.com/square/crossfilter/

*/

/*
    Utility functions, private to this module
*/

// fill an array or typedarray with a sequential range of numbers,
// starting with `start`
//
function fillRange(arr, start = 0) {
  for (let i = 0, len = arr.length; i < len; i++) {
    arr[i] = i + start;
  }
  return arr;
}

// Search for `value` in the sorted array `tarr`, in the range [first, last).
// Return the first (left most) index where tarr[index] >= value.
//
// In other words, return array index I where:
//  tarr[i] < value for all tarr[lo:I]
//  tarr[i] >= value for all tarr[I:last]
//
// Essentially the same thing as:
//    C++: lower_bound()
//    Python: bisect.bisect_left()
//
function lowerBound(valueArray, value, first, last) {
  // this is just a binary search
  while (first < last) {
    const middle = (first + last) >>> 1;
    if (valueArray[middle] < value) {
      first = middle + 1;
    } else {
      last = middle;
    }
  }
  return first;
}

// XXX: it is likely that there would be minimal performance hit from creating
// a factory version of lowerBound that takes an accessor (rather than having
// a special-cased version for lining the indirection).
//
// Benchmarking shows this manual inlining is up to 4X faster than an accessor.
// The real issue is how often we call it.
//
function lowerBoundIndirect(valueArray, indexArray, value, first, last) {
  // this is just a binary search
  while (first < last) {
    const middle = (first + last) >>> 1;
    if (valueArray[indexArray[middle]] < value) {
      first = middle + 1;
    } else {
      last = middle;
    }
  }
  return first;
}

// Search for `value in the sorted array `tarr`, in the range [first, last).
// Return the first value where tarr[index] > value.
//
// In other words, return array index I, where:
//  tarr[i] <= value for all tarr[lo:I]
//  tarr[i] > value for all tarr[I:last]
//
// Essentially the same thing as:
//    C++: upper_bound()
//    Python: bisect.bisect_right()
//
function upperBound(valueArray, value, first, last) {
  // this is just a binary search
  while (first < last) {
    const middle = (first + last) >>> 1;
    if (valueArray[middle] > value) {
      last = middle;
    } else {
      first = middle + 1;
    }
  }
  return first;
}

function upperBoundIndirect(valueArray, indexArray, value, first, last) {
  // this is just a binary search
  while (first < last) {
    const middle = (first + last) >>> 1;
    if (valueArray[indexArray[middle]] > value) {
      last = middle;
    } else {
      first = middle + 1;
    }
  }
  return first;
}

// Interval operations - very simple version of interval set relationship
// operators.   An interval is a multi-interval list of [min, max],
// where min and max are mandatory.   Constraints:
//    * min <= max
//    * Legal intervals:   [], [ [0, 1], ... ]
//    * all min and max values must be >= 0
//    * Not legal:   [ [] ]
//
// Code assumes intervals have a low cardinality; many operations are done
// with a brute force scan.
//
class PositiveIntervals {
  // Canonicalize - ensure that:
  //  1. no overlapping intervals
  //  2. sorted in order of interval min.
  //
  static canonicalize(A) {
    if (A.length <= 1) return A;
    let copy = A.slice();
    copy.sort((a, b) => a[0] - b[0]);
    const res = [];
    res.push(copy[0]);
    for (let i = 1, len = copy.length; i < len; i++) {
      if (copy[i][0] > res[res.length - 1][1]) {
        // non-overlapping, add to result
        res.push(copy[i]);
      } else if (copy[i][1] > res[res.length - 1][1]) {
        // merge this into previous
        res[res.length - 1][1] = copy[i][1];
      }
    }
    return res;
  }

  // Return interval with values belonging to both A and B.
  //
  static union(A, B) {
    return PositiveIntervals.canonicalize([...A, ...B]);
  }

  static _flatten(A, B) {
    let points = []; /* point, A, start */
    for (let a = 0; a < A.length; a++) {
      points.push([A[a][0], true, true]);
      points.push([A[a][1], true, false]);
    }
    for (let b = 0; b < B.length; b++) {
      points.push([B[b][0], false, true]);
      points.push([B[b][1], false, false]);
    }
    // Sort order: point, then start
    points.sort((a, b) => (a[0] !== b[0] ? a[0] - b[0] : a[2] ? 1 : -1));
    return points;
  }

  // A - B, ie, the interval with all values in A that are not in B.
  //
  static difference(A, B) {
    // Corner cases
    if (A.length === 0 || B.length === 0) {
      return PositiveIntervals.canonicalize(A);
    }

    A = PositiveIntervals.canonicalize(A);
    B = PositiveIntervals.canonicalize(B);

    const points = PositiveIntervals._flatten(A, B);
    const res = [];
    let aDepth = 0;
    let depth = 0;
    let intervalStart;
    let prevPoint;
    for (let i = 0; i < points.length; i++) {
      const p = points[i];
      const before = depth;
      const delta = p[2] ? 1 : -1;
      depth += delta;
      if (p[1]) aDepth += delta;

      if (i === points.length - 1 || p[0] !== points[i + 1][0]) {
        if (aDepth === 1 && depth === 1) {
          intervalStart = p[0];
        } else if (intervalStart !== undefined) {
          res.push([intervalStart, p[0]]);
          intervalStart = undefined;
        }
      }
      prevPoint = p[0];
    }
    // guaranteed to be in canonical form
    return res;
  }

  // Return interval with values belonging to A or B.
  //
  static intersection(A, B) {
    if (A.length === 0 || B.length === 0) {
      return [];
    }

    A = PositiveIntervals.canonicalize(A);
    B = PositiveIntervals.canonicalize(B);

    const points = PositiveIntervals._flatten(A, B);
    const res = [];
    let depth = 0;
    let intervalStart;
    for (let i = 0; i < points.length; i++) {
      const p = points[i];
      const before = depth;
      depth += p[2] ? 1 : -1;
      if (depth === 2) {
        intervalStart = p[0];
      } else if (intervalStart !== undefined) {
        res.push([intervalStart, p[0]]);
        intervalStart = undefined;
      }
    }
    // guaranteed to be in canonical form
    return res;
  }
}

// BitArray is a 2D bitarray with size [length, nBitWidth].
// Each bit is referred to as a `dimension`.  Dimensions may be
// dynamically allocated and deallocated.  The overall length
// of the BitArray is fixed at creation time (for simplicity).
//
// Organization of the bitarray is dimension-major.
//
// Primary operations on the BitArray are:
//    - set & clear dimension
//    - test dimension
//    - various performance or convenience test operations
//
// The underlying data structure uses TypedArrays for performance.
//
class BitArray {
  constructor(length) {
    // Initially allocate a 32 bit wide array.  allocDimension() will expand
    // as necessary.
    //
    // Int32Array is (counterintuitively) used to accomadate JS numeric casting
    // (to/from primitive number type).
    //

    // Fixed for the life of this object.
    this.length = length;

    // Bitarray width.  width is always greater than 32*dimensionCount.
    this.width = 1; // underlying number of 32 bit arrays
    this.dimensionCount = 0; // num allocated dimensions

    this.bitmask = new Int32Array(this.width); // dimension allocation mask
    this.bitarray = new Int32Array(this.width * this.length);
  }

  get selectionCount() {
    return this.countAllOnes();
  }

  countAllOnes() {
    let count = 0;
    for (let i = 0; i < this.width; i++) {
      const bitmask = this.bitmask[i];
      for (let j = i * this.length, len = j + this.length; j < len; j++) {
        if (this.bitarray[i * this.length + j] === bitmask) count++;
      }
    }
    return count;
  }

  // count trailing zeros
  static ctz(v) {
    let c = 32;
    v &= -v; // isolate lowest non-zero bit
    if (v) c--;
    if (v & 0x0000ffff) c -= 16;
    if (v & 0x00ff00ff) c -= 8;
    if (v & 0x0f0f0f0f) c -= 4;
    if (v & 0x33333333) c -= 2;
    if (v & 0x55555555) c -= 1;
    return c;
  }

  // find a free dimension.  Return undefined if none
  _findFreeDimension() {
    let dim;
    for (let col = 0; col < this.width; col++) {
      const bitmask = this.bitmask[col];
      const lowestZeroBit = ~this.bitmask[col] & -~this.bitmask[col];
      if (lowestZeroBit) {
        this.bitmask[col] |= lowestZeroBit;
        dim = 32 * col + BitArray.ctz(lowestZeroBit);
      }
    }
    return dim;
  }

  // allocate and return the dimension ID (bit position)
  allocDimension() {
    let dim = this._findFreeDimension();

    // if we did not find free dimension, expand the bitarray.
    if (dim === undefined) {
      this.width++;

      const biggerBitArray = new Int32Array(this.width * this.length);
      biggerBitArray.set(this.bitarray);
      this.bitarray = biggerBitArray;

      const biggerBitmask = new Int32Array(this.width);
      biggerBitmask.set(this.bitmask);
      this.bitmask = biggerBitmask;

      dim = this._findFreeDimension();
    }

    this.dimensionCount++;
    return dim;
  }

  freeDimension(dim) {
    // all selection tests assume unallocated dimensions are zero valued.
    this.deselectAll(dim);
    const col = dim >>> 5;
    this.bitmask[col] &= ~(1 << (dim % 32));
    this.dimensionCount--;
  }

  isSelected(index) {
    const width = this.width;
    const length = this.length;
    const bitarray = this.bitarray;

    for (let w = 0; w < width; w++) {
      const bitmask = this.bitmask[w];
      if (!bitmask || bitarray[w * length + index] !== bitmask) return false;
    }
    return true;
  }

  selectOne(dim, index) {
    const col = dim >>> 5;
    const before = this.bitarray[col * this.length + index];
    const after = before | (1 << (dim % 32));
    this.bitarray[col] = after;
  }

  deselectOne(dim, index) {
    const col = dim >>> 5;
    const before = this.bitarray[col * this.length + index];
    const after = before & ~(1 << (dim % 32));
    this.bitarray[col] = after;
  }

  selectAll(dim) {
    let col = dim >> 5;
    const bitmask = this.bitmask[col];
    const bitarray = this.bitarray;
    const one = 1 << (dim % 32);
    for (let i = col * this.length, len = i + this.length; i < len; i++) {
      bitarray[i] |= one;
    }
  }

  deselectAll(dim) {
    let col = dim >> 5;
    const bitmask = this.bitmask[col];
    const bitarray = this.bitarray;
    const zero = ~(1 << (dim % 32));
    for (let i = col * this.length, len = i + this.length; i < len; i++) {
      bitarray[i] &= zero;
    }
  }

  // indirect functions are used to map between sort and natural order
  selectIndirectFromRange(dim, indirect, range) {
    const col = dim >>> 5;
    const first = range[0];
    const last = range[1];
    const bitarray = this.bitarray;
    const one = 1 << (dim % 32);
    const offset = col * this.length;
    for (let i = first; i < last; i++) {
      bitarray[offset + indirect[i]] |= one;
    }
  }

  deselectIndirectFromRange(dim, indirect, range) {
    const col = dim >>> 5;
    const first = range[0];
    const last = range[1];
    const bitarray = this.bitarray;
    const zero = ~(1 << (dim % 32));
    const offset = col * this.length;
    for (let i = first; i < last; i++) {
      bitarray[offset + indirect[i]] &= zero;
    }
  }

  // Fill the array with selected|deselected value based upon the
  // current selection state.
  fillBySelection(result, selectedValue, deselectedValue) {
    // special case (width === 1) for performance
    if (this.width === 1) {
      const bitmask = this.bitmask[0];
      const bitarray = this.bitarray;
      for (let i = 0, len = this.length; i < len; i++) {
        result[i] = bitarray[i] === bitmask ? selectedValue : deselectedValue;
      }
    } else {
      for (let i = 0, len = this.length; i < len; i++) {
        result[i] = this.isSelected(i) ? selectedValue : deselectedValue;
      }
    }
    return result;
  }
}

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
    this.index = fillRange(new Uint32Array(this.crossfilter.data.length));
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
