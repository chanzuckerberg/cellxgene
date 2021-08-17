/*
Private dataframe support functions

TODO / XXX: for scalar/continuous data, this uses a naive method
of computing quantiles.  Would be good to switch from sort to 
partition at some point.
*/
import {
  AnyArray,
  GenericArrayConstructor,
} from "../../common/types/arraytypes";
import { ContinuousColumnSummary, CategoricalColumnSummary } from "./types";
import quantile from "../quantile";
import { sortArray } from "../typedCrossfilter/sort";

// [ 0, 0.01, 0.02, ..., 1.0]
const centileNames = new Array(101).fill(0).map((_v, idx) => idx / 100);

export function summarizeContinuous(col: AnyArray): ContinuousColumnSummary {
  let nan = 0;
  let pinf = 0;
  let ninf = 0;

  // -Inf < finite < Inf < NaN
  const sortedCol = sortArray(
    new (col.constructor as GenericArrayConstructor<typeof col>)(col)
  );

  // count non-finites, which are at each end of sorted data
  for (let i = sortedCol.length - 1; i >= 0; i -= 1) {
    if (!Number.isNaN(sortedCol[i])) {
      nan = sortedCol.length - i - 1;
      break;
    }
  }
  for (let i = 0, l = sortedCol.length; i < l; i += 1) {
    if (sortedCol[i] !== Number.NEGATIVE_INFINITY) {
      ninf = i;
      break;
    }
  }
  for (let i = sortedCol.length - nan - 1; i >= 0; i -= 1) {
    if (sortedCol[i] !== Number.POSITIVE_INFINITY) {
      pinf = sortedCol.length - i - nan - 1;
      break;
    }
  }

  // compute percentiles on finite data ONLY
  const sortedColFiniteOnly = sortedCol.slice(
    ninf,
    sortedCol.length - nan - pinf
  );
  const percentiles = quantile(centileNames, sortedColFiniteOnly, true);
  const min = percentiles[0];
  const max = percentiles[100];

  return {
    categorical: false,
    min,
    max,
    nan,
    pinf,
    ninf,
    percentiles,
  };
}

export function summarizeCategorical(col: AnyArray): CategoricalColumnSummary {
  const categoryCounts = new Map();
  if (col) {
    for (let r = 0, l = col.length; r < l; r += 1) {
      const val = col[r];
      let curCount = categoryCounts.get(val);
      if (curCount === undefined) curCount = 0;
      categoryCounts.set(val, curCount + 1);
    }
  }
  const sortedCategoryByCounts = new Map(
    [...categoryCounts.entries()].sort((a, b) => b[1] - a[1])
  );
  return {
    categorical: true,
    categories: [...sortedCategoryByCounts.keys()],
    categoryCounts: sortedCategoryByCounts,
    numCategories: sortedCategoryByCounts.size,
  };
}
