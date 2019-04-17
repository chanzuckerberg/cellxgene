/*
quantiles - calculate quantiles for the typed array.

Currently interpolates to 'higher' value.

Arguments:

	* q - array of quantiles to compute.
	* tarr - a typed array
*/

import { sort } from "./typedCrossfilter/sort";

export default function quantile(qs, tarr, sorted = false) {
	/*
	start with the naive (sort) implementation.  Later, use a faster partition
	*/
	const arr = sorted ? tarr : sort(new tarr.constructor(tarr)); // copy
	const len = arr.length;
	return qs.map(q => {
		if (q === 1) {
			return arr[len - 1];
		}
		return arr[Math.floor(q * len)];
	});
}
