/*
quantiles - calculate quantiles for the typed array.

Currently interpolates to 'lower' value.

Arguments:

	* quantArr - array of quantiles to compute, where values:  0 <= value <= 1.0
	* tarr - a typed array
	* sorted - option bool.  If false (default), will assume array is not sorted.
	 	If true, will assume it is sorted.

*/

import { sortArray } from "./typedCrossfilter/sort";

export default function quantile(
  quantArr: Array<number>,
  tarr:
    | Int8Array
    | Uint8Array
    | Int16Array
    | Uint16Array
    | Int32Array
    | Uint32Array
    | Uint8ClampedArray
    | Float32Array
    | Float64Array,
  sorted = false
): Array<unknown> {
  /*
	start with the naive (sort) implementation.  Later, use a faster partition
	*/
  const arr = sorted ? tarr : sortArray(tarr.constructor(tarr)); // copy
  const len = arr.length;
  return quantArr.map((q) => {
    if (q === 1) {
      return arr[len - 1];
    }
    return arr[Math.floor(q * len)];
  });
}
