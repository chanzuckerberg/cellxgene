import { sortIndex } from "./sort";
import { rangeFill as fillRange } from "../range";

/*
    Utility functions, private to this module.
*/

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
