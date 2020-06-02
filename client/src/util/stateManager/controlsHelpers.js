/*
Helper functions for the controls reducer
*/

import _ from "lodash";

import * as globals from "../../globals";
import { rangeFill as fillRange } from "../range";
import {
  userDefinedDimensionName,
  diffexpDimensionName,
} from "../nameCreators";

export function maxCategoryItems(config) {
  return (
    config.parameters?.["max-category-items"] ??
    globals.configDefaults.parameters["max-category-items"]
  );
}

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
  /* return top N by occurrences in the data, preserving original category order */
  const { categories } = colSchema;
  const counts = categories.map((cat) => summary.categoryCounts.get(cat) ?? 0);

  if (categories.length <= N) {
    return [categories, counts];
  }

  const sortIndex = fillRange(new Array(categories.length)).sort(
    (a, b) => counts[b] - counts[a]
  );
  const topNindices = new Set(sortIndex.slice(0, N));

  const _topNCategories = [];
  const topNCounts = [];
  for (let i = 0; i < categories.length; i += 1) {
    if (topNindices.has(i)) {
      _topNCategories.push(categories[i]);
      topNCounts.push(counts[i]);
    }
  }
  return [_topNCategories, topNCounts];
}

export function selectableCategoryNames(schema, maxCatItems, names) {
  /*
  return all obs annotation names that are categorical AND have a
  "reasonably" small number of categories AND are not the index column.

  If the initial name list not provided, use everything in the schema.
  */
  if (!schema) return [];
  const { index, columns } = schema.annotations.obs;

  return columns
    .filter((colSchema) => !names || names.indexOf(colSchema.name) !== -1)
    .filter((colSchema) => {
      const { type, name } = colSchema;
      const isSelectableType =
        type === "string" || type === "boolean" || type === "categorical";
      return isSelectableType && name !== index;
    })
    .map((v) => v.name);
}

export function createCategoricalSelection(world, names) {
  const N = globals.maxCategoricalOptionsToDisplay;
  const { obsAnnotations, schema } = world;

  const res = names.reduce((acc, name) => {
    const colSchema = schema.annotations.obsByName[name];
    const { writable: isUserAnno } = colSchema;

    /*
    Summarize the annotation data currently in world.  Must return categoryValues
    in sorted order, and must include all category values even if they are not
    actively used in the current world.
    */
    const summary = obsAnnotations.col(name).summarizeCategorical();
    const [categoryValues, categoryValueCounts] = topNCategories(
      colSchema,
      summary,
      N
    );
    const categoryValueIndices = new Map(categoryValues.map((v, i) => [v, i]));
    const numCategoryValues = categoryValueIndices.size;
    const categoryValueSelected = new Array(numCategoryValues).fill(true);
    const isTruncated = categoryValues.length < summary.numCategories;

    acc[name] = {
      categoryValues, // array: of natively typed category values
      categoryValueIndices, // map: category value (native type) -> category index
      categoryValueSelected, // array: t/f selection state
      numCategoryValues, // number: of values in the category
      isTruncated, // bool: true if list was truncated
      categoryValueCounts, // array: cardinality of each category,
      categorySelected: true, // bool - default state for entire category
      isUserAnno, // bool
    };
    return acc;
  }, {});
  return res;
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
