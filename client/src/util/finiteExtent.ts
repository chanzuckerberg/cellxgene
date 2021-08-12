/*
Return the [minimum, maximum] extent, of the given typed array, ignoring
non-finite values (ie, +Infinity, -Infinity).

If undefined or empty array, or array contains only non-finite numbers,
will return [undefined, undefined]
*/

import type { TypedArray } from "../common/types/entities";

function finiteExtent(
  tarr: TypedArray
): [number | undefined, number | undefined] {
  let min;
  let max;
  let i;

  for (i = 0; i < tarr.length; i += 1) {
    const val = tarr[i];
    if (Number.isFinite(val)) {
      min = val;
      max = val;
      i += 1;
      break;
    }
  }
  if (min && max)
    for (; i < tarr.length; i += 1) {
      const val = tarr[i];
      if (Number.isFinite(val)) {
        if (min > val) min = val;
        if (max < val) max = val;
      }
    }
  return [min, max];
}

export default finiteExtent;
