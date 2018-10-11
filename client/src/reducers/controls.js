// jshint esversion: 6

import _ from "lodash";
import { World, kvCache } from "../util/stateManager";
import { parseRGB } from "../util/parseRGB";
import Crossfilter from "../util/typedCrossfilter";
import * as globals from "../globals";
import {
  layoutDimensionName,
  obsAnnoDimensionName,
  userDefinedDimensionName,
  diffexpDimensionName,
  makeContinuousDimensionName
} from "../util/nameCreators";

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
    userDefinedGenes: [],
    diffexpGenes: [],

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
      const { userDefinedGenes, diffexpGenes } = state;
      /* first light - create world & other data-driven defaults */
      const { universe } = action;
      const world = World.createWorldFromEntireUniverse(universe);
      const colorName = new Array(universe.nObs).fill(globals.defaultCellColor);
      const colorRGB = _.map(colorName, c => parseRGB(c));
      const categoricalAsBooleansMap = createCategoricalAsBooleansMap(world);
      const crossfilter = Crossfilter(world.obsAnnotations);
      const dimensionMap = World.createObsDimensionMap(crossfilter, world);

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
        colorName,
        colorRGB,
        categoricalAsBooleansMap,
        crossfilter,
        dimensionMap,
        colorAccessor: null
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
      const colorName = new Array(world.nObs).fill(globals.defaultCellColor);
      const colorRGB = _.map(colorName, c => parseRGB(c));
      const categoricalAsBooleansMap = createCategoricalAsBooleansMap(world);
      const crossfilter = Crossfilter(world.obsAnnotations);
      const dimensionMap = World.createObsDimensionMap(crossfilter, world);

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
        colorName,
        colorRGB,
        categoricalAsBooleansMap,
        crossfilter,
        dimensionMap,
        colorAccessor: null
      };
    }
    case "expression load success": {
      console.log("expression load success", action);
      const { world, universe, crossfilter, dimensionMap } = state;
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

      _.forEach(action.expressionData, (val, key) => {
        dimensionMap[key] = World.createVarDimension(
          /* "__var__" + */
          world,
          worldVarDataCache,
          crossfilter,
          key
        );
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
    case "request differential expression success": {
      console.log("request diff success", action);
      const { world } = state;
      const _diffexpGenes = [];

      action.data.forEach(d => {
        _diffexpGenes.push(world.varAnnotations[d[0]].name);
      });

      return {
        ...state,
        diffexpGenes: _diffexpGenes
      };
    }
    case "clear differential expression": {
      const { world, universe, crossfilter, dimensionMap } = state;
      const _dimensionMap = dimensionMap;
      const universeVarDataCache = universe.varDataCache;
      const worldVarDataCache = world.varDataCache;

      _.forEach(action.diffExp, values => {
        const { name } = world.varAnnotations[values[0]].name;
        // clean up crossfilter dimensions
        const filterID = dimensionMap[diffexpDimensionName(name)];
        const dimension = crossfilter.filters.find(d => d.id === filterID);
        dimension.dim.dispose();

        // clean up dimensionsMap
        delete _dimensionMap[diffexpDimensionName(name)];
        // clean up the varDataCaches
        delete universeVarDataCache[name];
        delete worldVarDataCache[name];
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
      const { userDefinedGenes } = state;
      const newUserDefinedGenes = _.filter(
        userDefinedGenes,
        d => d !== action.data
      );
      return {
        ...state,
        userDefinedGenes: newUserDefinedGenes
      };
    }
    case "reset colorscale": {
      const { world } = state;
      const colorName = new Array(world.nObs).fill(globals.defaultCellColor);
      const colorRGB = _.map(colorName, c => parseRGB(c));
      return {
        ...state,
        colorName,
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
      state.dimensionMap[layoutDimensionName("X")].filterRange([
        action.brushCoords.northwest[0],
        action.brushCoords.southeast[0]
      ]);
      state.dimensionMap[layoutDimensionName("Y")].filterRange([
        action.brushCoords.southeast[1],
        action.brushCoords.northwest[1]
      ]);
      return {
        ...state,
        graphBrushSelection: action.brushCoords
      };
    }
    case "graph brush deselect": {
      state.dimensionMap[layoutDimensionName("X")].filterAll();
      state.dimensionMap[layoutDimensionName("Y")].filterAll();
      return {
        ...state,
        graphBrushSelection: null
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
      state.dimensionMap[obsAnnoDimensionName(action.metadataField)].filterEnum(
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
      state.dimensionMap[obsAnnoDimensionName(action.metadataField)].filterEnum(
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
      state.dimensionMap[
        obsAnnoDimensionName(action.metadataField)
      ].filterNone();
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
      state.dimensionMap[
        obsAnnoDimensionName(action.metadataField)
      ].filterAll();
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
