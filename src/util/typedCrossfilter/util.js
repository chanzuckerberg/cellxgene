"use strict";
// jshint esversion: 6

/*
    Utility functions, private to this module.
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

// Inlined performance optimization - used to indirect through a sort map.
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

// Inline performance optimization
//
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

module.exports = {
  fillRange,
  lowerBound,
  lowerBoundIndirect,
  upperBound,
  upperBoundIndirect
};
