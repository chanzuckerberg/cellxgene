// jshint esversion: 6

import _ from "lodash";

import { World, WorldUtil, ControlsHelper } from "../util/stateManager";
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

const Controls = (
  state = {
    // data loading flag
    loading: false,
    error: null,

    // configuration
    maxCategoryItems: globals.configDefaults.parameters["max-category-items"],

    // the whole big bang
    universe: null,
    fullUniverseCache: null,

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
    case "initial data load complete (universe exists)": {
      /* first light - create world & other data-driven defaults */
      const { universe } = action;
      const world = World.createWorldFromEntireUniverse(universe);
      const colorRGB = new Array(universe.nObs).fill(
        parseRGB(globals.defaultCellColor)
      );
      const categoricalSelectionState = ControlsHelper.createCategoricalSelectionState(
        state,
        world
      );
      const crossfilter = Crossfilter(world.obsAnnotations);
      const dimensionMap = World.createObsDimensionMap(crossfilter, world);
      WorldUtil.clearCaches();

      return {
        ...state,
        loading: false,
        error: null,
        universe,
        fullUniverseCache: { world, crossfilter, dimensionMap },
        world,
        colorRGB,
        categoricalSelectionState,
        crossfilter,
        dimensionMap,
        colorAccessor: null,
        resettingInterface: false
      };
    }
    case "reset World to eq Universe": {
      const {
        userDefinedGenes,
        diffexpGenes,
        universe,
        fullUniverseCache
      } = state;
      const { world, crossfilter } = fullUniverseCache;
      // reset all crossfilter dimensions
      _.forEach(fullUniverseCache.dimensionMap, dim => dim.filterAll());
      const colorRGB = new Array(universe.nObs).fill(
        parseRGB(globals.defaultCellColor)
      );
      const categoricalSelectionState = ControlsHelper.createCategoricalSelectionState(
        state,
        world
      );

      /* free dimensions not in cache (otherwise they leak) */
      _.forEach(state.dimensionMap, (dim, dimName) => {
        if (!fullUniverseCache.dimensionMap[dimName]) {
          dim.dispose();
        }
      });
      const dimensionMap = {
        ...fullUniverseCache.dimensionMap,
        ...ControlsHelper.createGenesDimMap(
          userDefinedGenes,
          diffexpGenes,
          world,
          crossfilter
        )
      };
      WorldUtil.clearCaches();

      return {
        ...state,
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
      const categoricalSelectionState = ControlsHelper.createCategoricalSelectionState(
        state,
        world
      );
      const crossfilter = Crossfilter(world.obsAnnotations);
      const dimensionMap = {
        ...World.createObsDimensionMap(crossfilter, world),
        ...ControlsHelper.createGenesDimMap(
          userDefinedGenes,
          diffexpGenes,
          world,
          crossfilter
        )
      };
      WorldUtil.clearCaches();

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
      let universeVarData = universe.varData;
      let worldVarData = world.varData;

      // Load new expression data into the varData dataframes, if
      // not already present.
      _.forEach(action.expressionData, (val, key) => {
        // If not already in universe.varData, save entire expression column
        if (!universeVarData.hasCol(key)) {
          universeVarData = universeVarData.withCol(key, val);
        }

        // If not already in world.varData, save sliced expression column
        if (!worldVarData.hasCol(key)) {
          // Slice if world !== universe, else just use whole column.
          // Use the obsAnnotation index as the cut key, as we keep
          // all world dataframes in sync.
          let worldValSlice = val;
          if (!World.worldEqUniverse(world, universe)) {
            worldValSlice = universeVarData
              .subset(world.obsAnnotations.rowIndex.keys(), [key], null)
              .icol(0)
              .asArray();
          }

          // Now build world's varData dataframe
          worldVarData = worldVarData.withCol(
            key,
            worldValSlice,
            world.obsAnnotations.rowIndex
          );
        }
      });

      // Prune size of varData "cache" if getting out of hand....
      const { userDefinedGenes, diffexpGenes } = state;
      const allTheGenesWeNeed = _.uniq(
        [].concat(
          userDefinedGenes,
          diffexpGenes,
          Object.keys(action.expressionData)
        )
      );
      universeVarData = ControlsHelper.pruneVarDataCache(
        universeVarData,
        allTheGenesWeNeed
      );
      worldVarData = ControlsHelper.pruneVarDataCache(
        worldVarData,
        allTheGenesWeNeed
      );

      return {
        ...state,
        universe: {
          ...universe,
          varData: universeVarData
        },
        world: {
          ...world,
          varData: worldVarData
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
      const _userDefinedGenes = userDefinedGenes.slice();
      const gene = action.data.genes[0];

      dimensionMap[
        userDefinedDimensionName(gene)
      ] = World.createVarDataDimension(world, crossfilter, gene);

      return {
        ...state,
        dimensionMap,
        userDefinedGenes: _userDefinedGenes,
        userDefinedGenesLoading: false
      };
    }
    case "request differential expression success": {
      const { world, crossfilter, dimensionMap } = state;
      const _diffexpGenes = [];

      action.data.forEach(d => {
        _diffexpGenes.push(world.varAnnotations.at(d[0], "name"));
      });

      _.forEach(_diffexpGenes, gene => {
        dimensionMap[diffexpDimensionName(gene)] = World.createVarDataDimension(
          world,
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
      const { world, dimensionMap } = state;
      const _dimensionMap = dimensionMap;
      _.forEach(action.diffExp, values => {
        const name = world.varAnnotations.at(values[0], "name");
        // clean up crossfilter dimensions
        const dimension = dimensionMap[diffexpDimensionName(name)];
        dimension.dispose();
        delete dimensionMap[diffexpDimensionName(name)];
      });
      return {
        ...state,
        dimensionMap: _dimensionMap,
        diffexpGenes: []
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
        ControlsHelper.selectedValuesForCategory(cat)
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
        ControlsHelper.selectedValuesForCategory(cat)
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
