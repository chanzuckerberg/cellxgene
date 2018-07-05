"use strict";
// jshint esversion: 6

// BitArray is a 2D bitarray with size [length, nBitWidth].
// Each bit is referred to as a `dimension`.  Dimensions may be
// dynamically allocated and deallocated.  The overall length
// of the BitArray is fixed at creation time (for simplicity).
//
// Organization of the bitarray is dimension-major.  As dimensions
// are added, the underlying store is grown 32 bits at a time.
// NOTE: currently does not deallocate / shrink.
//
// Primary operations on the BitArray are:
//    - set & clear dimension
//    - test dimension
//    - various performance or convenience operations to optimize bulk ops
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

  // Return the number of records that are selected, ie, have a one bit in
  // all allocated dimensions.
  //
  get selectionCount() {
    return this.countAllOnes();
  }

  // Count all records that have a 'one' bit in allocated dimensions.
  //
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

  // count trailing zeros - hard to do fast in JS!
  // https://en.wikipedia.org/wiki/Find_first_set#CTZ
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
  //
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

  // free a dimension for later use.  MUST deselect the dimension, as other
  // code assume the column will be zero valued.
  //
  freeDimension(dim) {
    // all selection tests assume unallocated dimensions are zero valued.
    this.deselectAll(dim);
    const col = dim >>> 5;
    this.bitmask[col] &= ~(1 << (dim % 32));
    this.dimensionCount--;
  }

  // return true if this index is selected in ALL dimensions.
  //
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

  // select index on dimension
  //
  selectOne(dim, index) {
    const col = dim >>> 5;
    const before = this.bitarray[col * this.length + index];
    const after = before | (1 << (dim % 32));
    this.bitarray[col * this.length + index] = after;
  }

  // deselect index on dimension
  //
  deselectOne(dim, index) {
    const col = dim >>> 5;
    const before = this.bitarray[col * this.length + index];
    const after = before & ~(1 << (dim % 32));
    this.bitarray[col * this.length + index] = after;
  }

  // select all indices on dimension.
  //
  selectAll(dim) {
    let col = dim >> 5;
    const bitmask = this.bitmask[col];
    const bitarray = this.bitarray;
    const one = 1 << (dim % 32);
    for (let i = col * this.length, len = i + this.length; i < len; i++) {
      bitarray[i] |= one;
    }
  }

  // deselect all indices on dimension
  //
  deselectAll(dim) {
    let col = dim >> 5;
    const bitmask = this.bitmask[col];
    const bitarray = this.bitarray;
    const zero = ~(1 << (dim % 32));
    for (let i = col * this.length, len = i + this.length; i < len; i++) {
      bitarray[i] &= zero;
    }
  }

  // select range of indices on a dimension, indirect through a sort map.
  // Indirect functions are used to map between sort and natural order.
  //
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

  // deselect range of indices on a dimension, indirect through a sort map.
  //
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
  //
  fillBySelection(result, selectedValue, deselectedValue) {
    // special case (width === 1) for performance
    if (this.width === 1) {
      const bitmask = this.bitmask[0];
      const bitarray = this.bitarray;
      for (let i = 0, len = this.length; i < len; i++) {
        result[i] =
          bitmask && bitarray[i] === bitmask ? selectedValue : deselectedValue;
      }
    } else {
      for (let i = 0, len = this.length; i < len; i++) {
        result[i] = this.isSelected(i) ? selectedValue : deselectedValue;
      }
    }
    return result;
  }
}

export default BitArray;
