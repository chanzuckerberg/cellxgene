import { sortIndex } from "./sort";
import { rangeFill as fillRange } from "../range";

/*
    Utility functions, private to this module.
*/

// slice out of one array into another, using an index array
//
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function sliceByIndex(src: any, index: any) {
  if (index === undefined || index === null) {
    return src;
  }
  const dst = new src.constructor(index.length);
  for (let i = 0; i < index.length; i += 1) {
    dst[i] = src[index[i]];
  }
  return dst;
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function makeSortIndex(src: any) {
  const index = fillRange(new Uint32Array(src.length));
  sortIndex(index, src);
  return index;
}
