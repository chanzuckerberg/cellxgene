// jshint esversion: 6

import { sortIndex } from "./sort";

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
  sortIndex(index, src);
  return index;
}
