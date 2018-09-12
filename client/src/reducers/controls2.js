// jshint esversion: 6

import _ from "lodash";
import { World } from "../util/stateManager";
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

    // all of the data + selection state
    world: null,
    categoricalAsBooleansMap: null, // should this be in world, w/ crossfilter?

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
      return Object.assign({}, state, {
        loading: true
      });
    }
    case "initial data load complete (universe exists)": {
      /* first light - create world & other data-driven defaults */
      const world = new World(action.universe, globals.defaultCellColor);
      const categoricalAsBooleansMap = createCategoricalAsBooleansMap(world);
      return Object.assign({}, state, {
        loading: false,
        error: null,
        world,
        categoricalAsBooleansMap,
        colorAccessor: null
      });
    }
    case "reset World to eq Universe": {
      /* Reset viewable world to be the entire Universe */
      const world = new World(state.world._universe, globals.defaultCellColor);
      const categoricalAsBooleansMap = createCategoricalAsBooleansMap(world);
      return Object.assign({}, state, {
        loading: false,
        error: null,
        world,
        categoricalAsBooleansMap,
        colorAccessor: null
      });
    }
    case "set World to current selection": {
      /* Set viewable world to be the currnetly selected data */
      const world = World.createFromCrossfilterSelection(
        state.world,
        globals.defaultCellColor
      );
      const categoricalAsBooleansMap = createCategoricalAsBooleansMap(world);
      return Object.assign({}, state, {
        world,
        categoricalAsBooleansMap,
        colorAccessor: null
      });
    }
    /* /api/v0.1/expression response */
    case "expression load success": {
      /* this is wrong and will be reworked when we use POJOs */
      state.world._universe.initFromExpression(action.data);
      return Object.assign({}, state);
    }
    case "expression load error":
    case "initial data load error": {
      return Object.assign({}, state, {
        loading: false,
        error: action.error
      });
    }

    /*******************************
             User Events
     *******************************/
    case "parallel coordinates axes have been drawn": {
      return Object.assign({}, state, {
        axesHaveBeenDrawn: true
      });
    }
    case "continuous selection using parallel coords brushing": {
      return Object.assign({}, state, {
        continuousSelection: action.data
      });
    }
    case "graph brush selection change": {
      state.world.obsDimensionMap.x.filterRange([
        action.brushCoords.northwest[0],
        action.brushCoords.southeast[0]
      ]);
      state.world.obsDimensionMap.y.filterRange([
        action.brushCoords.southeast[1],
        action.brushCoords.northwest[1]
      ]);
      return Object.assign({}, state, {
        graphBrushSelection: action.brushCoords
      });
    }
    case "graph brush deselect": {
      state.world.obsDimensionMap.x.filterAll();
      state.world.obsDimensionMap.y.filterAll();
      return Object.assign({}, state, {
        graphBrushSelection: null
      });
    }
    case "continuous metadata histogram brush": {
      // action.selection: metadata name being selected
      // action.range: filter range, or null if deselected
      if (!action.range) {
        state.world.obsDimensionMap[action.selection].filterAll();
      } else {
        state.world.obsDimensionMap[action.selection].filterRange(action.range);
      }
      return { ...state };
    }
    case "change opacity deselected cells in 2d graph background":
      return Object.assign({}, state, {
        opacityForDeselectedCells: action.data
      });

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
      state.world.obsDimensionMap[action.metadataField].filterEnum(
        _.filter(
          _.map(
            newCategoricalAsBooleansMap[action.metadataField],
            (val, key) => (val ? key : false)
          )
        )
      );
      return Object.assign({}, state, {
        categoricalAsBooleansMap: newCategoricalAsBooleansMap
      });
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
      state.world.obsDimensionMap[action.metadataField].filterEnum(
        _.filter(
          _.map(
            newCategoricalAsBooleansMap[action.metadataField],
            (val, key) => (val ? key : false)
          )
        )
      );
      return Object.assign({}, state, {
        categoricalAsBooleansMap: newCategoricalAsBooleansMap
      });
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
      state.world.obsDimensionMap[action.metadataField].filterNone();
      return Object.assign({}, state, {
        categoricalAsBooleansMap: newCategoricalAsBooleansMap
      });
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
      state.world.obsDimensionMap[action.metadataField].filterAll();
      return Object.assign({}, state, {
        categoricalAsBooleansMap: newCategoricalAsBooleansMap
      });
    }

    /*******************************
              Color Scale
    *******************************/
    case "color by categorical metadata":
    case "color by continuous metadata": {
      const world = state.world
        .clone()
        .setColors(action.colors.name, action.colors.rgb);
      return Object.assign({}, state, {
        colorAccessor: action.colorAccessor,
        colorScale: action.colorScale,
        world
      });
    }
    case "color by expression": {
      const world = state.world
        .clone()
        .setColors(action.colors.name, action.colors.rgb);
      return Object.assign({}, state, {
        colorAccessor: action.gene,
        colorScale: action.colorScale,
        world
      });
    }

    /*******************************
              Scatterplot
    *******************************/
    case "set scatterplot x":
      return Object.assign({}, state, {
        scatterplotXXaccessor: action.data
      });
    case "set scatterplot y":
      return Object.assign({}, state, {
        scatterplotYYaccessor: action.data
      });

    default:
      return state;
  }
};

export default Controls;
