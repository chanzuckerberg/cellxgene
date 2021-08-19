/*
Dataframe histogram
*/
import { isTypedArray } from "./util";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function _histogramContinuous(column: any, bins: any, min: any, max: any) {
  const valBins = new Array(bins).fill(0);
  if (!column) {
    return valBins;
  }
  const binWidth = (max - min) / bins;
  const colArray = column.asArray();
  for (let r = 0, len = colArray.length; r < len; r += 1) {
    const val = colArray[r];
    if (val <= max && val >= min) {
      // ensure test excludes NaN values
      const valBin = Math.min(Math.floor((val - min) / binWidth), bins - 1);
      valBins[valBin] += 1;
    }
  }
  return valBins;
}

function _histogramContinuousBy(
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  column: any,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  bins: any,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  min: any,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  max: any,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  by: any
) {
  const byMap = new Map();
  if (!column || !by) {
    return byMap;
  }
  const binWidth = (max - min) / bins;
  const byArray = by.asArray();
  const colArray = column.asArray();
  for (let r = 0, len = colArray.length; r < len; r += 1) {
    const byBin = byArray[r];
    let valBins = byMap.get(byBin);
    if (valBins === undefined) {
      valBins = new Array(bins).fill(0);
      byMap.set(byBin, valBins);
    }
    const val = colArray[r];
    if (val <= max && val >= min) {
      // ensure test excludes NaN values
      const valBin = Math.min(Math.floor((val - min) / binWidth), bins - 1);
      valBins[valBin] += 1;
    }
  }
  return byMap;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function _histogramCategorical(column: any) {
  const valMap = new Map();
  if (!column) {
    return valMap;
  }
  const colArray = column.asArray();
  for (let r = 0, len = colArray.length; r < len; r += 1) {
    const valBin = colArray[r];
    let curCount = valMap.get(valBin);
    if (curCount === undefined) {
      curCount = 0;
    }
    valMap.set(valBin, curCount + 1);
  }
  return valMap;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function _histogramCategoricalBy(column: any, by: any) {
  const byMap = new Map();
  if (!column || !by) {
    return byMap;
  }
  const byArray = by.asArray();
  const colArray = column.asArray();
  for (let r = 0, len = colArray.length; r < len; r += 1) {
    const byBin = byArray[r];
    let valMap = byMap.get(byBin);
    if (valMap === undefined) {
      valMap = new Map();
      byMap.set(byBin, valMap);
    }
    const valBin = colArray[r];
    let curCount = valMap.get(valBin);
    if (curCount === undefined) {
      curCount = 0;
    }
    valMap.set(valBin, curCount + 1);
  }
  return byMap;
}

/*
Count category occupancy.  Optional group-by category.
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function histogramCategorical(column: any, by: any) {
  if (by && isTypedArray(by)) {
    throw new Error("Group by column must be categorical");
  }
  return by
    ? _histogramCategoricalBy(column, by)
    : _histogramCategorical(column);
}

/*
Memoization hash for histogramCategorical()
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function hashCategorical(column: any, by: any) {
  if (by) {
    return `${column.__id}:${by.__id}`;
  }
  return `${column.__id}:`;
}

/*
Bin counts for continuous/scalar values, with optional group-by category.
Values outside domain are ignored.
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function histogramContinuous(
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  column: any,
  bins = 40,
  domain = [0, 1],
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  by: any
) {
  if (by && isTypedArray(by)) {
    throw new Error("Group by column must be categorical");
  }
  const [min, max] = domain;
  return by
    ? _histogramContinuousBy(column, bins, min, max, by)
    : _histogramContinuous(column, bins, min, max);
}

/*
Memoization hash for histogramContinuous
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function hashContinuous(
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  column: any,
  bins = "",
  domain = [0, 0],
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  by: any
) {
  const [min, max] = domain;
  if (by) {
    return `${column.__id}:${bins}:${min}:${max}:${by.__id}`;
  }
  return `${column.__id}::${bins}:${min}:${max}`;
}
