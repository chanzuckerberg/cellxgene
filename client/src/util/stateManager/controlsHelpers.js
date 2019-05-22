/*
Helper functions for the controls reducer
*/

import _ from "lodash";

import * as globals from "../../globals";
import { rangeFill as fillRange } from "../range";
import {
  userDefinedDimensionName,
  diffexpDimensionName
} from "../nameCreators";

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
function topNCategories(summary) {
  const counts = _.map(summary.categories, cat =>
    summary.categoryCounts.get(cat)
  );
  const sortIndex = fillRange(new Array(summary.numCategories)).sort(
    (a, b) => counts[b] - counts[a]
  );
  const sortedCategories = _.map(sortIndex, i => summary.categories[i]);
  const sortedCounts = _.map(sortIndex, i => counts[i]);
  const N = globals.maxCategoricalOptionsToDisplay;

  if (sortedCategories.length < N) {
    return [sortedCategories, sortedCounts];
  }
  return [sortedCategories.slice(0, N), sortedCounts.slice(0, N)];
}

export function createCategoricalSelection(maxCategoryItems, world) {
  const res = {};
  _.forEach(world.obsAnnotations.colIndex.keys(), key => {
    const summary = world.obsAnnotations.col(key).summarize();
    if (summary.categories) {
      const isColorField = key.includes("color") || key.includes("Color");
      const isSelectableCategory =
        !isColorField &&
        key !== "name" &&
        summary.categories.length < maxCategoryItems;
      if (isSelectableCategory) {
        const [categoryValues, categoryValueCounts] = topNCategories(summary);
        const categoryValueIndices = new Map(
          categoryValues.map((v, i) => [v, i])
        );
        const numCategoryValues = categoryValueIndices.size;
        const categoryValueSelected = new Array(numCategoryValues).fill(true);
        const isTruncated = categoryValues.length < summary.numCategories;
        res[key] = {
          categoryValues, // array: of natively typed category values
          categoryValueIndices, // map: category value (native type) -> category index
          categoryValueSelected, // array: t/f selection state
          numCategoryValues, // number: of values in the category
          isTruncated, // bool: true if list was truncated
          categoryValueCounts, // array: cardinality of each category,
          categorySelected: true // bool - default state for entire category
        };
      }
    }
  });
  return res;
}

/*
given a categoricalSelection, return the list of all category values
where selection state is true (ie, they are selected).
*/
export function selectedValuesForCategory(categorySelectionState, dfColumn) {
  const {
    categorySelected,
    categoryValueSelected,
    categoryValueIndices
  } = categorySelectionState;
  let selectedValues;
  if (categorySelected) {
    selectedValues = new Set(dfColumn.summarize().categories);
  } else {
    selectedValues = new Set();
  }
  categoryValueIndices.forEach((catIndex, catValue) => {
    if (!categoryValueSelected[catIndex]) {
      selectedValues.delete(catValue);
    } else {
      selectedValues.add(catValue);
    }
  });
  return [...selectedValues.values()];
}

/*
build a crossfilter dimensions for all gene expression related dimensions.
*/
export function createGeneDimensions(
  userDefinedGenes,
  diffexpGenes,
  world,
  crossfilter
) {
  crossfilter = userDefinedGenes.reduce(
    (xflt, gene) =>
      xflt.addDimension(
        userDefinedDimensionName(gene),
        "scalar",
        world.varData.col(gene).asArray(),
        Float32Array
      ),
    crossfilter
  );
  crossfilter = diffexpGenes.reduce(
    (xflt, gene) =>
      xflt.addDimension(
        diffexpDimensionName(gene),
        "scalar",
        world.varData.col(gene).asArray(),
        Float32Array
      ),
    crossfilter
  );
  return crossfilter;
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
  VarDataCacheLowWatermark - this cofig value sets the minimum cache size,
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
  const all = colIndex.keys();
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
