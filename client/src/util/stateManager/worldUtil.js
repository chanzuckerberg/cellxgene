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

Parameters are:
- dim1: dimension 1 name/label
- dim2: dimension 2 name/label
- df: dataframe containing dim1 and dim2 on the column axis

*/
function _countCategoryValues2D(dim1, dim2, df) {
  const dimMap = new Map();
  const col1 = df.col(dim1) ? df.col(dim1).asArray() : null;
  const col2 = df.col(dim2) ? df.col(dim2).asArray() : null;
  if (!col1 || !col2) {
    return dimMap;
  }

  for (let r = 0, l = df.length; r < l; r += 1) {
    const val1 = col1[r];
    const val2 = col2[r];
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

let __worldUtilMemoId__ = 0;
function _memoizedId(x) {
  if (!x.__worldUtilMemoId__) {
    __worldUtilMemoId__ += 1;
    x.__worldUtilMemoId__ = __worldUtilMemoId__;
  }
  return x.__worldUtilMemoId__;
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
Clear any cached data within WorldUtil caches, eg, memoized functions
*/
export function clearCaches() {
  countCategoryValues2D.cache.clear();
}
