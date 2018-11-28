/* eslint-disable import/prefer-default-export */
import _ from "lodash";

/*
Various utility functions operating on World/Universe
*/

/*
Count unique category values, binning first by dim1 then by dim2
Return:

Map {
   dim1_val1: Map {
      dim2_val1: number,
      dim2_val2: number,
      ...
   },
   ...
}

*/
function _countCategoryValues2D(dim1, dim2, rows) {
  const dimMap = new Map();
  for (let r = 0; r < rows.length; r += 1) {
    const row = rows[r];
    const val1 = row[dim1];
    const val2 = row[dim2];
    let d2Map = dimMap.get(val1);
    if (d2Map === undefined) {
      d2Map = new Map();
      dimMap.set(val1, d2Map);
    }
    let curCount = d2Map.get(val2);
    if (curCount === undefined) {
      curCount = 0;
    }
    d2Map.set(val2, curCount + 1);
  }
  return dimMap;
}

let __worldOpsMemoId__ = 0;
function _memoizedId(x) {
  if (!x.__worldOpsMemoId__) {
    __worldOpsMemoId__ += 1;
    x.__worldOpsMemoId__ = __worldOpsMemoId__;
  }
  return x.__worldOpsMemoId__;
}
function _countCategoryValues2DResolver(...args) {
  const id = args[0] + args[1] + _memoizedId(args[2]);
  return id;
}

export const countCategoryValues2D = _.memoize(
  _countCategoryValues2D,
  _countCategoryValues2DResolver
);

/*
Clear any cached data within WorldOps caches, eg, memoized functions
*/
export function clearCaches() {
  countCategoryValues2D.cache.clear();
}
