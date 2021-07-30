/* eslint-disable no-bitwise -- crossfilter relies on bitwise ops */

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
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  bitarray: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  bitmask: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  dimensionCount: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  length: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  width: any;

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  constructor(length: any) {
    // Initially allocate a 32 bit wide array.  allocDimension() will expand
    // as necessary.
    //
    // Int32Array is (counterintuitively) used to accomadate JS numeric casting
    // (to/from primitive number type).
    //

    // Fixed for the life of this object.
    this.length = length;

    // Bitarray width.  width is always greater than dimensionCount/32.
    this.width = 1; // underlying number of 32 bit arrays
    this.dimensionCount = 0; // num allocated dimensions

    this.bitmask = new Int32Array(this.width); // dimension allocation mask
    this.bitarray = new Int32Array(this.width * this.length);
    Object.seal(this);
  }

  // Return the number of records that are selected, ie, have a one bit in
  // all allocated dimensions.
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  selectionCount() {
    return this.countAllOnes();
  }

  // Count all records that have a 'one' bit in allocated dimensions.
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  countAllOnes() {
    let count = 0;
    const { bitarray, length, width } = this;
    if (width === 1) {
      // special case, width === 1, for performance
      const bitmask = this.bitmask[0];
      for (let l = 0; l < length; l += 1) {
        if (bitarray[l] === bitmask) {
          count += 1;
        }
      }
    } else {
      const { bitmask } = this;
      for (let l = 0; l < length; l += 1) {
        let dimensionsSet = 0;
        for (let w = 0; w < width; w += 1) {
          if (bitarray[w * length + l] === bitmask[w]) {
            dimensionsSet += 1;
          }
        }
        if (dimensionsSet === width) {
          count += 1;
        }
      }
    }
    return count;
  }

  // count trailing zeros - hard to do fast in JS!
  // https://en.wikipedia.org/wiki/Find_first_set#CTZ
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  static ctz(av: any) {
    let c = 32;
    let v = av;
    v &= -v; // isolate lowest non-zero bit
    if (v) c -= 1;
    if (v & 0x0000ffff) c -= 16;
    if (v & 0x00ff00ff) c -= 8;
    if (v & 0x0f0f0f0f) c -= 4;
    if (v & 0x33333333) c -= 2;
    if (v & 0x55555555) c -= 1;
    return c;
  }

  // find a free dimension.  Return undefined if none
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  _findFreeDimension() {
    let dim;
    for (let col = 0; col < this.width; col += 1) {
      const lowestZeroBit = ~this.bitmask[col] & -~this.bitmask[col];
      if (lowestZeroBit) {
        this.bitmask[col] |= lowestZeroBit;
        dim = 32 * col + BitArray.ctz(lowestZeroBit);
        break;
      }
    }
    return dim;
  }

  // allocate and return the dimension ID (bit position)
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  allocDimension() {
    let dim = this._findFreeDimension();

    // if we did not find free dimension, expand the bitarray.
    if (dim === undefined) {
      this.width += 1;

      const biggerBitArray = new Int32Array(this.width * this.length);
      biggerBitArray.set(this.bitarray);
      this.bitarray = biggerBitArray;

      const biggerBitmask = new Int32Array(this.width);
      biggerBitmask.set(this.bitmask);
      this.bitmask = biggerBitmask;

      dim = this._findFreeDimension();
    }

    this.dimensionCount += 1;
    return dim;
  }

  // free a dimension for later use.  MUST deselect the dimension, as other
  // code assume the column will be zero valued.
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  freeDimension(dim: any) {
    // all selection tests assume unallocated dimensions are zero valued.
    this.deselectAll(dim);
    const col = dim >>> 5;
    this.bitmask[col] &= ~(1 << dim % 32);
    this.dimensionCount -= 1;
  }

  // return true if this index is selected in ALL dimensions.
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  isSelected(index: any) {
    const { width, length, bitarray } = this;

    for (let w = 0; w < width; w += 1) {
      const bitmask = this.bitmask[w];
      if (!bitmask || bitarray[w * length + index] !== bitmask) return false;
    }
    return true;
  }

  // return true if this index is selected in ALL dimensions IGNORING dim
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  isSelectedIgnoringDim(index: any, dim: any) {
    const ignoreOffset = dim >>> 5;
    const ignoreMask = ~(1 << dim % 32);

    const { width, length, bitarray } = this;

    for (let w = 0; w < width; w += 1) {
      const bitmask = this.bitmask[w];
      if (w === ignoreOffset) {
        if (
          bitmask &&
          (bitarray[w * length + index] & ignoreMask) !== (bitmask & ignoreMask)
        ) {
          return false;
        }
      } else if (bitmask && bitarray[w * length + index] !== bitmask) {
        return false;
      }
    }
    return true;
  }

  // select index on dimension
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  selectOne(dim: any, index: any) {
    const col = dim >>> 5;
    const before = this.bitarray[col * this.length + index];
    const after = before | (1 << dim % 32);
    this.bitarray[col * this.length + index] = after;
  }

  // deselect index on dimension
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  deselectOne(dim: any, index: any) {
    const col = dim >>> 5;
    const before = this.bitarray[col * this.length + index];
    const after = before & ~(1 << dim % 32);
    this.bitarray[col * this.length + index] = after;
  }

  // select all indices on dimension.
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  selectAll(dim: any) {
    const col = dim >> 5;
    const one = 1 << dim % 32;
    for (let i = col * this.length, len = i + this.length; i < len; i += 1) {
      this.bitarray[i] |= one;
    }
  }

  // deselect all indices on dimension
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  deselectAll(dim: any) {
    const col = dim >> 5;
    const zero = ~(1 << dim % 32);
    for (let i = col * this.length, len = i + this.length; i < len; i += 1) {
      this.bitarray[i] &= zero;
    }
  }

  // select range of indices on a dimension
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  selectFromRange(dim: any, range: any) {
    const col = dim >>> 5;
    const first = range[0];
    const last = range[1];
    const one = 1 << dim % 32;
    const offset = col * this.length;
    for (let i = first; i < last; i += 1) {
      this.bitarray[offset + i] |= one;
    }
  }

  // select range of indices on a dimension, indirect through a sort map.
  // Indirect functions are used to map between sort and natural order.
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  selectIndirectFromRange(dim: any, indirect: any, range: any) {
    const col = dim >>> 5;
    const first = range[0];
    const last = range[1];
    const one = 1 << dim % 32;
    const offset = col * this.length;
    for (let i = first; i < last; i += 1) {
      this.bitarray[offset + indirect[i]] |= one;
    }
  }

  // deselect range of indices on a dimension
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  deselectFromRange(dim: any, range: any) {
    const col = dim >>> 5;
    const first = range[0];
    const last = range[1];
    const zero = ~(1 << dim % 32);
    const offset = col * this.length;
    for (let i = first; i < last; i += 1) {
      this.bitarray[offset + i] &= zero;
    }
  }

  // deselect range of indices on a dimension, indirect through a sort map.
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  deselectIndirectFromRange(dim: any, indirect: any, range: any) {
    const col = dim >>> 5;
    const first = range[0];
    const last = range[1];
    const zero = ~(1 << dim % 32);
    const offset = col * this.length;
    for (let i = first; i < last; i += 1) {
      this.bitarray[offset + indirect[i]] &= zero;
    }
  }

  // Fill the array with selected|deselected value based upon the
  // current selection state.
  //
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  fillBySelection(result: any, selectedValue: any, deselectedValue: any) {
    // special case (width === 1) for performance
    if (this.width === 1) {
      const { bitmask, bitarray } = this;
      const mask = bitmask[0];
      if (!mask) {
        result.fill(deselectedValue);
      } else {
        for (let i = 0, len = this.length; i < len; i += 1) {
          result[i] = bitarray[i] === mask ? selectedValue : deselectedValue;
        }
      }
    } else {
      for (let i = 0, len = this.length; i < len; i += 1) {
        result[i] = this.isSelected(i) ? selectedValue : deselectedValue;
      }
    }
    return result;
  }
}

export default BitArray;
/* eslint-enable no-bitwise -- enable */
