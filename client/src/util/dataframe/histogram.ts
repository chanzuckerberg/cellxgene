/*
Dataframe histogram
*/
import { NumberArray } from "../../common/types/arraytypes";
import {
  DataframeColumn,
  ContinuousHistogram,
  ContinuousHistogramBy,
  CategoricalHistogram,
  CategoricalHistogramBy,
} from "./types";

export function histogramContinuous(
  column: DataframeColumn,
  bins: number,
  domain: [number, number]
): ContinuousHistogram {
  const valBins = new Array(bins).fill(0);
  if (!column) {
    return valBins;
  }
  const [min, max] = domain;
  const binWidth = (max - min) / bins;
  const colArray: NumberArray = column.asArray() as NumberArray;
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

export function histogramContinuousBy(
  column: DataframeColumn,
  bins: number,
  domain: [number, number],
  by: DataframeColumn
): ContinuousHistogramBy {
  const byMap = new Map();
  if (!column || !by) {
    return byMap;
  }
  const [min, max] = domain;
  const binWidth = (max - min) / bins;
  const byArray = by.asArray();
  const colArray = column.asArray() as NumberArray;
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

export function histogramCategorical(
  column: DataframeColumn
): CategoricalHistogram {
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

export function histogramCategoricalBy(
  column: DataframeColumn,
  by: DataframeColumn
): CategoricalHistogramBy {
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
Memoization hash for histogramCategorical()
*/
export function hashCategorical(
  column: DataframeColumn,
  by?: DataframeColumn
): string {
  if (by) {
    return `${column.__id}:${by.__id}`;
  }
  return `${column.__id}:`;
}

/*
Memoization hash for histogramContinuous
*/
export function hashContinuous(
  column: DataframeColumn,
  bins = "",
  domain: [number, number] = [0, 0],
  by?: DataframeColumn
): string {
  const [min, max] = domain;
  if (by) {
    return `${column.__id}:${bins}:${min}:${max}:${by.__id}`;
  }
  return `${column.__id}::${bins}:${min}:${max}`;
}
