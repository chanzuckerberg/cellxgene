// jshint esversion: 6

import _ from "lodash";
import { World } from "../util/stateManager";

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
  console.log(action);
  switch (action.type) {
    /*****************************************************
          Initialization and World/Universe management
    ******************************************************/
    case "initial data load complete (universe exists)":
    case "reset World to eq Universe": {
      /* Reset viewable world to be the entire Universe */
      /* first light - create world & other data-driven defaults */
      const world = new World(action.universe, "rgb(0,0,0,1)");
      const categoricalAsBooleansMap = createCategoricalAsBooleansMap(world);
      return Object.assign({}, state, {
        world,
        categoricalAsBooleansMap
      });
    }
    case "set World to current selection": {
      /* Set viewable world to be the currnetly selected data */
      const world = World.createFromCrossfilterSelection(state.world);
      const categoricalAsBooleansMap = createCategoricalAsBooleansMap(world);
      return Object.assign({}, state, {
        world,
        categoricalAsBooleansMap
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
