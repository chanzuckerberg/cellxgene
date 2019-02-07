// jshint esversion: 6

import _ from "lodash";
import { polygonContains } from "d3";

import { World, kvCache, WorldUtil } from "../util/stateManager";
import parseRGB from "../util/parseRGB";
import Crossfilter from "../util/typedCrossfilter";
import * as globals from "../globals";
import {
  layoutDimensionName,
  obsAnnoDimensionName,
  userDefinedDimensionName,
  diffexpDimensionName,
  makeContinuousDimensionName
} from "../util/nameCreators";
import { fillRange } from "../util/typedCrossfilter/util";

/*
Selection state for categoricals are tracked in an Object that
has two main components for each category:
1. mapping of option value to an index
2. array of bool selection state by index
Remember that option values can be ANY js type, except undefined/null.

  {
    _category_name_1: {
      // map of option value to index
      categoryIndices: Map([
        catval1: index,
        ...
      ])

      // index->selection true/false state
      categorySelected: [ true/false, true/false, ... ]

      // number of options
      numCategories: number,

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

function createCategoricalSelectionState(state, world) {
  const res = {};
  _.forEach(world.summary.obs, (value, key) => {
    if (value.categories) {
      const isColorField = key.includes("color") || key.includes("Color");
      const isSelectableCategory =
        !isColorField &&
        key !== "name" &&
        value.categories.length < state.maxCategoryItems;
      if (isSelectableCategory) {
        const [categoryValues, categoryCounts] = topNCategories(value);
        const categoryIndices = new Map(categoryValues.map((v, i) => [v, i]));
        const numCategories = categoryIndices.size;
        const categorySelected = new Array(numCategories).fill(true);
        const isTruncated = categoryValues.length < value.numCategories;
        res[key] = {
          categoryValues, // array: of natively typed category values
          categoryIndices, // map: category value (native type) -> category index
          categorySelected, // array: t/f selection state
          numCategories, // number: of categories
          isTruncated, // bool: true if list was truncated
          categoryCounts // array: cardinality of each category
        };
      }
    }
  });
  return res;
}

/*
given a categoricalSelectionState, return the list of all category values
where selection state is true (ie, they are selected).
*/
function selectedValuesForCategory(categorySelectionState) {
  const selectedValues = _([...categorySelectionState.categoryIndices])
    .filter(tuple => categorySelectionState.categorySelected[tuple[1]])
    .map(tuple => tuple[0])
    .value();
  return selectedValues;
}

const Controls = (
  state = {
    // data loading flag
    loading: false,
    error: null,

    // configuration
    maxCategoryItems: globals.configDefaults.parameters["max-category-items"],

    // the whole big bang
    universe: null,

    // all of the data + selection state
    world: null,
    colorRGB: null,
    categoricalSelectionState: null,
    crossfilter: null,
    dimensionMap: null,
    userDefinedGenes: [],
    userDefinedGenesLoading: false,
    diffexpGenes: [],

    colorAccessor: null,
    colorScale: null,
    resettingInterface: false,

    opacityForDeselectedCells: 0.2,
    graphBrushSelection: null,
    continuousSelection: null,
    scatterplotXXaccessor: null, // just easier to read
    scatterplotYYaccessor: null,
    axesHaveBeenDrawn: false,
    graphRenderCounter: 0 /* integer as <Component key={graphRenderCounter} - a change in key forces a remount */,
    __storedStateForCelllist1__: null /* will need procedural control of brush ie., brush.extent https://bl.ocks.org/micahstubbs/3cda05ca68cba260cb81 */,
    __storedStateForCelllist2__: null
  },
  action
) => {
  /*
  For now, log anything looking like an error to the console.
  */
  if (action.error || /error/i.test(action.type)) {
    console.error(action.error);
  }

  switch (action.type) {
    /*****************************************************
          Initialization, World/Universe management
          and data loading.
    ******************************************************/
    case "configuration load complete": {
      // there are a couple of configuration items we need to retain
      return {
        ...state,
        maxCategoryItems: _.get(
          state.config,
          "parameters.max-category-items",
          globals.configDefaults.parameters["max-category-items"]
        )
      };
    }
    case "initial data load start": {
      return { ...state, loading: true };
    }
    case "initial data load complete (universe exists)":
    case "reset World to eq Universe": {
      const { userDefinedGenes, diffexpGenes } = state;
      /* first light - create world & other data-driven defaults */
      const { universe } = action;
      const world = World.createWorldFromEntireUniverse(universe);
      const colorRGB = new Array(universe.nObs).fill(
        parseRGB(globals.defaultCellColor)
      );
      const categoricalSelectionState = createCategoricalSelectionState(
        state,
        world
      );
      const crossfilter = Crossfilter(world.obsAnnotations);
      const dimensionMap = World.createObsDimensionMap(crossfilter, world);
      WorldUtil.clearCaches();

      const worldVarDataCache = world.varDataCache;

      // dimensionMap = {
      //     layout_X: dim-for-X,
      //     obsAnno_name: dim for an annotation,
      //     varData_userDefined_genename: dim for user defined expression,
      //     varData_diffexp_genename: dim for diff-exp added gene expression
      // }
      /* var dimensions */
      if (userDefinedGenes.length > 0) {
        /*
          verbose & slightly confusing that we also access this as an object
          in controls rather than an array, should be abstracted into
          util ie., createDimensionsFromBothListsOfGenes(userGenes, diffExp)
        */
        _.forEach(userDefinedGenes, gene => {
          dimensionMap[
            userDefinedDimensionName(gene)
          ] = World.createVarDimension(
            /* "__var__" + */
            world,
            worldVarDataCache,
            crossfilter,
            gene
          );
        });
      }

      if (diffexpGenes.length > 0) {
        _.forEach(diffexpGenes, gene => {
          dimensionMap[diffexpDimensionName(gene)] = World.createVarDimension(
            /* "__var__" + */
            world,
            worldVarDataCache,
            crossfilter,
            gene
          );
        });
      }

      return {
        ...state,
        loading: false,
        error: null,
        universe,
        world,
        colorRGB,
        categoricalSelectionState,
        crossfilter,
        dimensionMap,
        colorAccessor: null,
        resettingInterface: false
      };
    }
    case "set World to current selection": {
      const { userDefinedGenes, diffexpGenes } = state;

      /* Set viewable world to be the currently selected data */
      const world = World.createWorldFromCurrentSelection(
        action.universe,
        action.world,
        action.crossfilter
      );
      const colorRGB = new Array(world.nObs).fill(
        parseRGB(globals.defaultCellColor)
      );
      const categoricalSelectionState = createCategoricalSelectionState(
        state,
        world
      );
      const crossfilter = Crossfilter(world.obsAnnotations);
      const dimensionMap = World.createObsDimensionMap(crossfilter, world);
      WorldUtil.clearCaches();

      const worldVarDataCache = world.varDataCache;
      /* var dimensions */

      if (userDefinedGenes.length > 0) {
        /*
          verbose & slightly confusing that we also access this as an object
          in controls rather than an array, should be abstracted into
          util ie., createDimensionsFromBothListsOfGenes(userGenes, diffExp)
        */
        _.forEach(userDefinedGenes, gene => {
          dimensionMap[
            userDefinedDimensionName(gene)
          ] = World.createVarDimension(
            /* "__var__" + */
            world,
            worldVarDataCache,
            crossfilter,
            gene
          );
        });
      }

      if (diffexpGenes.length > 0) {
        _.forEach(diffexpGenes, gene => {
          dimensionMap[diffexpDimensionName(gene)] = World.createVarDimension(
            /* "__var__" + */
            world,
            worldVarDataCache,
            crossfilter,
            gene
          );
        });
      }

      return {
        ...state,
        loading: false,
        error: null,
        world,
        colorRGB,
        categoricalSelectionState,
        crossfilter,
        dimensionMap,
        colorAccessor: null
      };
    }
    case "expression load success": {
      const { world, universe } = state;
      let universeVarDataCache = universe.varDataCache;
      let worldVarDataCache = world.varDataCache;
      _.forEach(action.expressionData, (val, key) => {
        universeVarDataCache = kvCache.set(universeVarDataCache, key, val);
        if (kvCache.get(worldVarDataCache, key) === undefined) {
          worldVarDataCache = kvCache.set(
            worldVarDataCache,
            key,
            World.subsetVarData(world, universe, val)
          );
        }
      });

      return {
        ...state,
        universe: {
          ...universe,
          varDataCache: universeVarDataCache
        },
        world: {
          ...world,
          varDataCache: worldVarDataCache
        }
      };
    }
    case "request user defined gene started": {
      return {
        ...state,
        userDefinedGenesLoading: true
      };
    }
    case "request user defined gene error": {
      return {
        ...state,
        userDefinedGenesLoading: false
      };
    }
    case "request user defined gene success": {
      const { world, crossfilter, dimensionMap, userDefinedGenes } = state;
      const worldVarDataCache = world.varDataCache;
      const _userDefinedGenes = userDefinedGenes.slice();
      const gene = action.data.genes[0];

      dimensionMap[userDefinedDimensionName(gene)] = World.createVarDimension(
        /* "__var__" + */
        world,
        worldVarDataCache,
        crossfilter,
        gene
      );

      return {
        ...state,
        dimensionMap,
        userDefinedGenes: _userDefinedGenes,
        userDefinedGenesLoading: false
      };
    }
    case "request differential expression success": {
      const { world, crossfilter, dimensionMap } = state;
      const worldVarDataCache = world.varDataCache;
      const _diffexpGenes = [];

      action.data.forEach(d => {
        _diffexpGenes.push(world.varAnnotations[d[0]].name);
      });

      _.forEach(_diffexpGenes, gene => {
        dimensionMap[diffexpDimensionName(gene)] = World.createVarDimension(
          /* "__var__" + */
          world,
          worldVarDataCache,
          crossfilter,
          gene
        );
      });

      return {
        ...state,
        dimensionMap,
        diffexpGenes: _diffexpGenes
      };
    }
    case "clear differential expression": {
      const { world, universe, dimensionMap } = state;
      const _dimensionMap = dimensionMap;
      const universeVarDataCache = universe.varDataCache;
      const worldVarDataCache = world.varDataCache;

      _.forEach(action.diffExp, values => {
        const { name } = world.varAnnotations[values[0]];
        // clean up crossfilter dimensions
        const dimension = dimensionMap[diffexpDimensionName(name)];
        dimension.dispose();
        delete dimensionMap[diffexpDimensionName(name)];
      });
      return {
        ...state,
        dimensionMap: _dimensionMap,
        diffexpGenes: [],
        universe: {
          ...universe,
          varDataCache: universeVarDataCache
        },
        world: {
          ...world,
          varDataCache: worldVarDataCache
        }
      };
    }
    case "user defined gene": {
      /*
        this could also live in expression success with a conditional,
        but that handles diffexp also
      */
      const newUserDefinedGenes = state.userDefinedGenes.slice();
      newUserDefinedGenes.push(action.data);
      return {
        ...state,
        userDefinedGenes: newUserDefinedGenes
      };
    }
    case "clear user defined gene": {
      const { userDefinedGenes, dimensionMap } = state;
      const newUserDefinedGenes = _.filter(
        userDefinedGenes,
        d => d !== action.data
      );

      const dimension = dimensionMap[userDefinedDimensionName(action.data)];
      dimension.dispose();
      delete dimensionMap[userDefinedDimensionName(action.data)];

      return {
        ...state,
        dimensionMap,
        userDefinedGenes: newUserDefinedGenes
      };
    }
    case "clear all user defined genes": {
      const { userDefinedGenes, dimensionMap } = state;

      _.forEach(userDefinedGenes, gene => {
        const dimension = dimensionMap[userDefinedDimensionName(gene)];
        dimension.dispose();
        delete dimensionMap[userDefinedDimensionName(gene)];
      });

      return {
        ...state,
        dimensionMap,
        userDefinedGenes: []
      };
    }
    case "reset colorscale": {
      const { world } = state;
      const colorRGB = new Array(world.nObs).fill(
        parseRGB(globals.defaultCellColor)
      );
      return {
        ...state,
        colorRGB,
        colorAccessor: null
      };
    }
    case "expression load error":
    case "initial data load error": {
      return {
        ...state,
        loading: false,
        error: action.error
      };
    }

    /*******************************
             User Events
     *******************************/
    case "graph brush selection change": {
      state.dimensionMap[layoutDimensionName("XY")].filterWithinRect(
        action.brushCoords.northwest,
        action.brushCoords.southeast
      );
      return {
        ...state,
        graphBrushSelection: action.brushCoords
      };
    }
    case "lasso deselect":
    case "graph brush deselect": {
      state.dimensionMap[layoutDimensionName("XY")].filterAll();
      return {
        ...state,
        graphBrushSelection: null
      };
    }
    case "lasso selection": {
      const { polygon } = action;
      const dXY = state.dimensionMap[layoutDimensionName("XY")];
      if (polygon.length < 3) {
        // single point or a line is not a polygon, and is therefore a deselect
        dXY.filterAll();
      } else {
        dXY.filterWithinPolygon(polygon);
      }
      return {
        ...state
      };
    }
    case "continuous metadata histogram brush": {
      const name = makeContinuousDimensionName(
        action.continuousNamespace,
        action.selection
      );

      // action.selection: metadata name being selected
      // action.range: filter range, or null if deselected
      if (!action.range) {
        state.dimensionMap[name].filterAll();
      } else {
        state.dimensionMap[name].filterRange(action.range);
      }
      return { ...state };
    }
    case "change opacity deselected cells in 2d graph background":
      return {
        ...state,
        opacityForDeselectedCells: action.data
      };
    case "increment graph render counter": {
      const c = state.graphRenderCounter + 1;
      return {
        ...state,
        graphRenderCounter: c
      };
    }
    case "interface reset started": {
      return {
        ...state,
        resettingInterface: true
      };
    }
    /*******************************
          Categorical metadata
    *******************************/
    case "categorical metadata filter select": {
      const newCategorySelected = Array.from(
        state.categoricalSelectionState[action.metadataField].categorySelected
      );
      newCategorySelected[action.categoryIndex] = true;
      const newCategoricalSelectionState = {
        ...state.categoricalSelectionState,
        [action.metadataField]: {
          ...state.categoricalSelectionState[action.metadataField],
          categorySelected: newCategorySelected
        }
      };

      // update the filter to match all selected options
      const cat = newCategoricalSelectionState[action.metadataField];
      state.dimensionMap[obsAnnoDimensionName(action.metadataField)].filterEnum(
        selectedValuesForCategory(cat)
      );

      return {
        ...state,
        categoricalSelectionState: newCategoricalSelectionState
      };
    }
    case "categorical metadata filter deselect": {
      const newCategorySelected = Array.from(
        state.categoricalSelectionState[action.metadataField].categorySelected
      );
      newCategorySelected[action.categoryIndex] = false;
      const newCategoricalSelectionState = {
        ...state.categoricalSelectionState,
        [action.metadataField]: {
          ...state.categoricalSelectionState[action.metadataField],
          categorySelected: newCategorySelected
        }
      };

      // update the filter to match all selected options
      const cat = newCategoricalSelectionState[action.metadataField];
      state.dimensionMap[obsAnnoDimensionName(action.metadataField)].filterEnum(
        selectedValuesForCategory(cat)
      );

      return {
        ...state,
        categoricalSelectionState: newCategoricalSelectionState
      };
    }
    case "categorical metadata filter none of these": {
      const newCategoricalSelectionState = {
        ...state.categoricalSelectionState,
        [action.metadataField]: {
          ...state.categoricalSelectionState[action.metadataField],
          categorySelected: Array.from(
            state.categoricalSelectionState[action.metadataField]
              .categorySelected
          ).fill(false)
        }
      };
      state.dimensionMap[
        obsAnnoDimensionName(action.metadataField)
      ].filterNone();
      return {
        ...state,
        categoricalSelectionState: newCategoricalSelectionState
      };
    }
    case "categorical metadata filter all of these": {
      const newCategoricalSelectionState = {
        ...state.categoricalSelectionState,
        [action.metadataField]: {
          ...state.categoricalSelectionState[action.metadataField],
          categorySelected: Array.from(
            state.categoricalSelectionState[action.metadataField]
              .categorySelected
          ).fill(true)
        }
      };
      state.dimensionMap[
        obsAnnoDimensionName(action.metadataField)
      ].filterAll();
      return {
        ...state,
        categoricalSelectionState: newCategoricalSelectionState
      };
    }

    /*******************************
              Color Scale
    *******************************/
    case "color by categorical metadata":
    case "color by continuous metadata": {
      return {
        ...state,
        colorRGB: action.colors.rgb,
        colorAccessor: action.colorAccessor,
        colorScale: action.colorScale
      };
    }
    case "color by expression": {
      return {
        ...state,
        colorRGB: action.colors.rgb,
        colorAccessor: action.gene,
        colorScale: action.colorScale
      };
    }

    /*******************************
              Scatterplot
    *******************************/
    case "set scatterplot x":
      return {
        ...state,
        scatterplotXXaccessor: action.data
      };
    case "set scatterplot y":
      return {
        ...state,
        scatterplotYYaccessor: action.data
      };
    case "clear scatterplot":
      return {
        ...state,
        scatterplotXXaccessor: null,
        scatterplotYYaccessor: null
      };

    default:
      return state;
  }
};

export default Controls;
