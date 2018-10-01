// jshint esversion: 6

import _ from "lodash";
import { World, kvCache } from "../util/stateManager";
import { parseRGB } from "../util/parseRGB";
import Crossfilter from "../util/typedCrossfilter";
import * as globals from "../globals";

function createCategoricalAsBooleansMap(world) {
  const res = {};
  _.each(world.summary.obs, (value, key) => {
    if (value.options) {
      const optionsAsBooleans = {};
      _.each(value.options, (_value, _key) => {
        optionsAsBooleans[_key] = true;
      });
      res[key] = optionsAsBooleans;
    }
  });
  return res;
}

const Controls = (
  state = {
    // data loading flag
    loading: false,
    error: null,

    universe: null,

    // all of the data + selection state
    world: null,
    colorName: null,
    colorRGB: null,
    categoricalAsBooleansMap: null,
    crossfilter: null,
    dimensionMap: null,
    dimensionMapVar: {},

    colorAccessor: null,
    colorScale: null,

    opacityForDeselectedCells: 0.2,
    graphBrushSelection: null,
    continuousSelection: null,
    scatterplotXXaccessor: null, // just easier to read
    scatterplotYYaccessor: null,
    axesHaveBeenDrawn: false,
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
    case "initial data load start": {
      return { ...state, loading: true };
    }
    case "initial data load complete (universe exists)":
    case "reset World to eq Universe": {
      /* first light - create world & other data-driven defaults */
      const { universe } = action;
      const world = World.createWorldFromEntireUniverse(universe);
      const colorName = new Array(universe.nObs).fill(globals.defaultCellColor);
      const colorRGB = _.map(colorName, c => parseRGB(c));
      const categoricalAsBooleansMap = createCategoricalAsBooleansMap(world);
      const crossfilter = Crossfilter(world.obsAnnotations);
      const dimensionMap = World.createObsDimensionMap(crossfilter, world);
      return {
        ...state,
        loading: false,
        error: null,
        universe,
        world,
        colorName,
        colorRGB,
        categoricalAsBooleansMap,
        crossfilter,
        dimensionMap,
        colorAccessor: null
      };
    }
    case "set World to current selection": {
      /* Set viewable world to be the currently selected data */
      const world = World.createWorldFromCurrentSelection(
        action.universe,
        action.world,
        action.crossfilter
      );
      const colorName = new Array(world.nObs).fill(globals.defaultCellColor);
      const colorRGB = _.map(colorName, c => parseRGB(c));
      const categoricalAsBooleansMap = createCategoricalAsBooleansMap(world);
      const crossfilter = Crossfilter(world.obsAnnotations);
      const dimensionMap = World.createObsDimensionMap(crossfilter, world);
      return {
        ...state,
        loading: false,
        error: null,
        world,
        colorName,
        colorRGB,
        categoricalAsBooleansMap,
        crossfilter,
        dimensionMap,
        colorAccessor: null
      };
    }
    case "expression load success": {
      const {
        world,
        universe,
        crossfilter,
        dimensionMap,
        dimensionMapVar
      } = state;
      let universeVarDataCache = universe.varDataCache;
      let worldVarDataCache = world.varDataCache;
      let _dimensionMapVar = dimensionMapVar;
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

      _.forEach(action.expressionData, (val, key) => {
        console.log("got one ", key);
        dimensionMap[key] = World.createVarDimensionsMap(
          world,
          worldVarDataCache,
          crossfilter,
          key
        );
      });

      return {
        ...state,
        dimensionMapVar,
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
      state.dimensionMap.x.filterRange([
        action.brushCoords.northwest[0],
        action.brushCoords.southeast[0]
      ]);
      state.dimensionMap.y.filterRange([
        action.brushCoords.southeast[1],
        action.brushCoords.northwest[1]
      ]);
      return {
        ...state,
        graphBrushSelection: action.brushCoords
      };
    }
    case "graph brush deselect": {
      state.dimensionMap.x.filterAll();
      state.dimensionMap.y.filterAll();
      return {
        ...state,
        graphBrushSelection: null
      };
    }
    case "continuous metadata histogram brush": {
      // action.selection: metadata name being selected
      // action.range: filter range, or null if deselected
      if (!action.range) {
        state.dimensionMap[action.selection].filterAll();
      } else {
        state.dimensionMap[action.selection].filterRange(action.range);
      }
      return { ...state };
    }
    case "change opacity deselected cells in 2d graph background":
      return {
        ...state,
        opacityForDeselectedCells: action.data
      };

    /*******************************
          Categorical metadata
    *******************************/
    case "categorical metadata filter select": {
      const newCategoricalAsBooleansMap = {
        ...state.categoricalAsBooleansMap,
        [action.metadataField]: {
          ...state.categoricalAsBooleansMap[action.metadataField],
          [action.value]: true
        }
      };
      // update the filter for the one category that changed state
      state.dimensionMap[action.metadataField].filterEnum(
        _.filter(
          _.map(
            newCategoricalAsBooleansMap[action.metadataField],
            (val, key) => (val ? key : false)
          )
        )
      );
      return {
        ...state,
        categoricalAsBooleansMap: newCategoricalAsBooleansMap
      };
    }
    case "categorical metadata filter deselect": {
      const newCategoricalAsBooleansMap = {
        ...state.categoricalAsBooleansMap,
        [action.metadataField]: {
          ...state.categoricalAsBooleansMap[action.metadataField],
          [action.value]: false
        }
      };
      // update the filter for the one category that changed state
      state.dimensionMap[action.metadataField].filterEnum(
        _.filter(
          _.map(
            newCategoricalAsBooleansMap[action.metadataField],
            (val, key) => (val ? key : false)
          )
        )
      );
      return {
        ...state,
        categoricalAsBooleansMap: newCategoricalAsBooleansMap
      };
    }
    case "categorical metadata filter none of these": {
      const newCategoricalAsBooleansMap = {
        ...state.categoricalAsBooleansMap
      };
      _.forEach(
        newCategoricalAsBooleansMap[action.metadataField],
        (v, k, c) => {
          c[k] = false;
        }
      );
      state.dimensionMap[action.metadataField].filterNone();
      return {
        ...state,
        categoricalAsBooleansMap: newCategoricalAsBooleansMap
      };
    }
    case "categorical metadata filter all of these": {
      const newCategoricalAsBooleansMap = {
        ...state.categoricalAsBooleansMap
      };
      _.forEach(
        newCategoricalAsBooleansMap[action.metadataField],
        (v, k, c) => {
          c[k] = true;
        }
      );
      state.dimensionMap[action.metadataField].filterAll();
      return {
        ...state,
        categoricalAsBooleansMap: newCategoricalAsBooleansMap
      };
    }

    /*******************************
              Color Scale
    *******************************/
    case "color by categorical metadata":
    case "color by continuous metadata": {
      return {
        ...state,
        colorName: action.colors.name,
        colorRGB: action.colors.rgb,
        colorAccessor: action.colorAccessor,
        colorScale: action.colorScale
      };
    }
    case "color by expression": {
      return {
        ...state,
        colorName: action.colors.name,
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

    default:
      return state;
  }
};

export default Controls;
