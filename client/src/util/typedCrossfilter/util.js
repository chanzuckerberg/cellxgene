// jshint esversion: 6
/* eslint no-bitwise: "off" */

/*
    Utility functions, private to this module.
*/

// fill an array or typedarray with a sequential range of numbers,
// starting with `start`
//
export function fillRange(arr, start = 0) {
  const larr = arr;
  for (let i = 0, len = larr.length; i < len; i += 1) {
    larr[i] = i + start;
  }
  return larr;
}

// slice out of one array into another, using an index array
//
export function sliceByIndex(src, index) {
  if (index === undefined || index === null) {
    return src;
  }
  const dst = new src.constructor(index.length);
  for (let i = 0; i < index.length; i += 1) {
    dst[i] = src[index[i]];
  }
  return dst;
}

export function makeSortIndex(src) {
  const index = fillRange(new Uint32Array(src.length));
  index.sort((a, b) => src[a] - src[b]);
  return index;
}

// Search for `value` in the sorted array `arr`, in the range [first, last).
// Return the first (left most) index where arr[index] >= value.
//
// In other words, return array index I where:
//  arr[i] < value for all tarr[lo:I]
//  arr[i] >= value for all tarr[I:last]
//
// The same semantics/behavior as:
//    C++: lower_bound()
//    Python: bisect.bisect_left()
//
// XXX: it is likely that there would be minimal performance hit from creating
// a factory version of lowerBound that takes an accessor (rather than having
// a special-cased version for lining the indirection).
//
export function lowerBound(valueArray, value, first, last) {
  let lfirst = first;
  let llast = last;
  // this is just a binary search
  while (lfirst < llast) {
    const middle = (lfirst + llast) >>> 1;
    if (valueArray[middle] < value) {
      lfirst = middle + 1;
    } else {
      llast = middle;
    }
  }
  return lfirst;
}

// Inlined performance optimization - used to indirect through a sort map.
//
export function lowerBoundIndirect(valueArray, indexArray, value, first, last) {
  let lfirst = first;
  let llast = last;
  // this is just a binary search
  while (lfirst < llast) {
    const middle = (lfirst + llast) >>> 1;
    if (valueArray[indexArray[middle]] < value) {
      lfirst = middle + 1;
    } else {
      llast = middle;
    }
  }
  return lfirst;
}

// Search for `value in the sorted array `arr`, in the range [first, last).
// Return the first value where arr[index] > value.
//
// In other words, return array index I, where:
//  arr[i] <= value for all tarr[lo:I]
//  arr[i] > value for all tarr[I:last]
//
// The same semantics/behavior as:
//    C++: upper_bound()
//    Python: bisect.bisect_right()
//
export function upperBound(valueArray, value, first, last) {
  let lfirst = first;
  let llast = last;
  // this is just a binary search
  while (lfirst < llast) {
    const middle = (lfirst + llast) >>> 1;
    if (valueArray[middle] > value) {
      llast = middle;
    } else {
      lfirst = middle + 1;
    }
  }
  return lfirst;
}

// Inline performance optimization
//
export function upperBoundIndirect(valueArray, indexArray, value, first, last) {
  let lfirst = first;
  let llast = last;
  // this is just a binary search
  while (lfirst < llast) {
    const middle = (lfirst + llast) >>> 1;
    if (valueArray[indexArray[middle]] > value) {
      llast = middle;
    } else {
      lfirst = middle + 1;
    }
  }
  return lfirst;
}
