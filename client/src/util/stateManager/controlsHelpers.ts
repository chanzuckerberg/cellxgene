/*
Helper functions for the controls reducer
*/

import difference from "lodash.difference";

import fromEntries from "../fromEntries";
import { isCategoricalAnnotation } from "./annotationsHelpers";

/*
Selection state for categoricals are tracked in an Object that
has two main components for each category:
1. mapping of option value to an index
2. array of bool selection state by index
Remember that option values can be ANY js type, except undefined/null.

  {
    _category_name_1: {
      // map of option value to index
      categoryValueIndices: Map([
        catval1: index,
        ...
      ])

      // index->selection true/false state
      categoryValueSelected: [ true/false, true/false, ... ]

      // number of options
      numCategoryValues: number,
    }
  }
*/

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function isSelectableCategoryName(schema: any, name: any) {
  const { index } = schema.annotations.obs;
  const colSchema = schema.annotations.obsByName[name];
  return (
    name &&
    name !== index &&
    (isCategoricalAnnotation(schema, name) || colSchema.writable)
  );
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function selectableCategoryNames(schema: any, names: any) {
  /*
  return all obs annotation names that are categorical AND have a
  "reasonably" small number of categories AND are not the index column.

  If the initial name list not provided, use everything in the schema.
  */
  if (!schema) return [];
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  if (!names) names = schema.annotations.obs.columns.map((c: any) => c.name);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  return names.filter((name: any) => isSelectableCategoryName(schema, name));
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function createCategorySummaryFromDfCol(dfCol: any, colSchema: any) {
  const { writable: isUserAnno } = colSchema;

  /*
  Summarize the annotation data currently in dataframe column.  Must return
  categoryValues in sorted order, and must include all category values even
  if they are not actively used in the current annoMatrix view.
  */
  const summary = dfCol.summarizeCategorical();
  const { categories: allCategoryValues } = colSchema;
  const categoryValues = allCategoryValues;
  const categoryValueCounts = allCategoryValues.map(
    (cat: any) => summary.categoryCounts.get(cat) ?? 0
  );
  const categoryValueIndices = new Map(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    categoryValues.map((v: any, i: any) => [v, i])
  );
  const numCategoryValues = categoryValueIndices.size;

  return {
    allCategoryValues, // array: of natively typed category values (all of them)
    categoryValues, // array: of natively typed category values (top N only)
    categoryValueIndices, // map: category value (native type) -> category index (top N only)
    numCategoryValues, // number: of values in the category (top N)
    categoryValueCounts, // array: cardinality of each category, (top N)
    isUserAnno, // bool
  };
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function createCategoricalSelection(names: any) {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  return fromEntries(names.map((name: any) => [name, new Map()]));
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function pruneVarDataCache(varData: any, needed: any) {
  /*
  Remove any unneeded columns from the varData dataframe.  Will only
  prune / remove if the total column count exceeds VarDataCacheLowWatermark

  Note: this code leverages the fact that dataframe offsets indicate
  the order in which the columns were added.  This crudely provides
  LRU semantics, so we can delete "older" columns first.
  */

  /*
  VarDataCacheLowWatermark - this config value sets the minimum cache size,
  in columns, below which we don't throw away data.

  The value should be high enough so we are caching the maximum which will
  "typically" be used in the UI (currently: 10 for diffexp, and N for user-
  specified genes), and low enough to account for memory use (any single
  column size is 4 bytes * numObs, so a column can be multi-megabyte in common
  use cases).
  */
  const VarDataCacheLowWatermark = 32;

  const numOverWatermark = varData.dims[1] - VarDataCacheLowWatermark;
  if (numOverWatermark <= 0) return varData;

  const { colIndex } = varData;
  const all = colIndex.labels();
  const unused = difference(all, needed);
  if (unused.length > 0) {
    // sort by offset in the dataframe - ie, psuedo-LRU
    unused.sort((a, b) => colIndex.getOffset(a) - colIndex.getOffset(b));
    const numToDrop =
      unused.length < numOverWatermark ? unused.length : numOverWatermark;
    for (let i = 0; i < numToDrop; i += 1) {
      varData = varData.dropCol(unused[i]);
    }
  }
  return varData;
}
