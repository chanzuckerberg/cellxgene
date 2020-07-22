/*
Helper functions for the controls reducer
*/

import _ from "lodash";

import * as globals from "../../globals";
import { rangeFill as fillRange } from "../range";
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

      // isTruncated - true if the options for selection has
      // been truncated (ie, was too large to implement)
    }
  }
*/
function topNCategories(colSchema, summary, N) {
  /* return top N categories by occurrences in the data */
  const { categories: allCategories } = colSchema;
  const counts = allCategories.map(
    (cat) => summary.categoryCounts.get(cat) ?? 0
  );

  if (allCategories.length <= N) {
    return [allCategories, allCategories, counts];
  }

  const sortIndex = fillRange(new Array(allCategories.length)).sort(
    (a, b) => counts[b] - counts[a]
  );
  const topNindices = new Set(sortIndex.slice(0, N));

  const _topNCategories = [];
  const topNCounts = [];
  for (let i = 0; i < allCategories.length; i += 1) {
    if (topNindices.has(i)) {
      _topNCategories.push(allCategories[i]);
      topNCounts.push(counts[i]);
    }
  }
  return [allCategories, _topNCategories, topNCounts];
}

export function isSelectableCategoryName(schema, name) {
  const { index } = schema.annotations.obs;
  const colSchema = schema.annotations.obsByName[name];
  return (
    name &&
    name !== index &&
    (isCategoricalAnnotation(schema, name) || colSchema.writable)
  );
}

export function selectableCategoryNames(schema, names) {
  /*
  return all obs annotation names that are categorical AND have a
  "reasonably" small number of categories AND are not the index column.

  If the initial name list not provided, use everything in the schema.
  */
  if (!schema) return [];
  if (!names) names = schema.annotations.obs.columns.map((c) => c.name);
  return names.filter((name) => isSelectableCategoryName(schema, name));
}

export function createCategorySummaryFromDfCol(dfCol, colSchema) {
  const N = globals.maxCategoricalOptionsToDisplay;
  const { writable: isUserAnno } = colSchema;

  /*
  Summarize the annotation data currently in dataframe column.  Must return
  categoryValues in sorted order, and must include all category values even
  if they are not actively used in the current annoMatrix view.
  */
  const summary = dfCol.summarizeCategorical();
  const [
    allCategoryValues,
    categoryValues,
    categoryValueCounts,
  ] = topNCategories(colSchema, summary, N);
  const categoryValueIndices = new Map(categoryValues.map((v, i) => [v, i]));
  const numCategoryValues = categoryValueIndices.size;
  const isTruncated = categoryValues.length < summary.numCategories;

  return {
    allCategoryValues, // array: of natively typed category values (all of them)
    categoryValues, // array: of natively typed category values (top N only)
    categoryValueIndices, // map: category value (native type) -> category index (top N only)
    numCategoryValues, // number: of values in the category (top N)
    isTruncated, // bool: true if list was truncated (ie, if topN != all)
    categoryValueCounts, // array: cardinality of each category, (top N)
    isUserAnno, // bool
  };
}

export function createCategoricalSelection(names) {
  return fromEntries(names.map((name) => [name, new Map()]));
}

export function pruneVarDataCache(varData, needed) {
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
  const unused = _.difference(all, needed);
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

export function subsetAndResetGeneLists(state) {
  const { userDefinedGenes, diffexpGenes } = state;
  const newUserDefinedGenes = _.uniq(
    [].concat(userDefinedGenes, diffexpGenes)
  ).slice(0, globals.maxGenes);
  const newDiffExpGenes = [];
  return [newUserDefinedGenes, newDiffExpGenes];
}
