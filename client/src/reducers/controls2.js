// jshint esversion: 6

import _ from "lodash";
import { World, kvCache } from "../util/stateManager";
import { parseRGB } from "../util/parseRGB";
import crossfilter from "../util/typedCrossfilter";
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
    loading: false, // XXX: separate reducer?
    error: null,

    universe: null, // XXX: separate reducer?

    // all of the data + selection state
    world: null,
    colorName: null,
    colorRGB: null,
    categoricalAsBooleansMap: null, // should this be in world, w/ crossfilter?
    obsCrossfilter: null, // separate reducer?
    obsDimensionMap: null, // separate reducer?

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
      const obsCrossfilter = crossfilter(world.obsAnnotations);
      const obsDimensionMap = World.createObsDimensionMap(
        obsCrossfilter,
        world
      );
      return {
        ...state,
        loading: false,
        error: null,
        universe,
        world,
        colorName,
        colorRGB,
        categoricalAsBooleansMap,
        obsCrossfilter,
        obsDimensionMap,
        colorAccessor: null
      };
    }
    case "set World to current selection": {
      /* Set viewable world to be the currently selected data */
      const world = World.createWorldFromCurrentSelection(
        action.universe,
        action.world,
        action.obsCrossfilter
      );
      const colorName = new Array(world.nObs).fill(globals.defaultCellColor);
      const colorRGB = _.map(colorName, c => parseRGB(c));
      const categoricalAsBooleansMap = createCategoricalAsBooleansMap(world);
      const obsCrossfilter = crossfilter(world.obsAnnotations);
      const obsDimensionMap = World.createObsDimensionMap(
        obsCrossfilter,
        world
      );
      return {
        ...state,
        loading: false,
        error: null,
        world,
        colorName,
        colorRGB,
        categoricalAsBooleansMap,
        obsCrossfilter,
        obsDimensionMap,
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
    case "parallel coordinates axes have been drawn": {
      return {
        ...state,
        axesHaveBeenDrawn: true
      };
    }
    case "continuous selection using parallel coords brushing": {
      return {
        ...state,
        continuousSelection: action.data
      };
    }
    case "graph brush selection change": {
      state.obsDimensionMap.x.filterRange([
        action.brushCoords.northwest[0],
        action.brushCoords.southeast[0]
      ]);
      state.obsDimensionMap.y.filterRange([
        action.brushCoords.southeast[1],
        action.brushCoords.northwest[1]
      ]);
      return {
        ...state,
        graphBrushSelection: action.brushCoords
      };
    }
    case "graph brush deselect": {
      state.obsDimensionMap.x.filterAll();
      state.obsDimensionMap.y.filterAll();
      return {
        ...state,
        graphBrushSelection: null
      };
    }
    case "continuous metadata histogram brush": {
      // action.selection: metadata name being selected
      // action.range: filter range, or null if deselected
      if (!action.range) {
        state.obsDimensionMap[action.selection].filterAll();
      } else {
        state.obsDimensionMap[action.selection].filterRange(action.range);
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
      state.obsDimensionMap[action.metadataField].filterEnum(
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
      state.obsDimensionMap[action.metadataField].filterEnum(
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
      state.obsDimensionMap[action.metadataField].filterNone();
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
      state.obsDimensionMap[action.metadataField].filterAll();
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
