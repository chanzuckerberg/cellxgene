// jshint esversion: 6
import _ from "lodash";
import { parseRGB } from "../util/parseRGB";

const Controls = (
  state = {
    _ranges: null /* this comes from initialize, this is universe */,
    allGeneNames: null,
    allCellsOnClient: null /* this comes from cells endpoint, this is world */,
    currentCellSelection: null /* this comes from user actions, all draw components use this, it is created by middleware */,
    graphMap: null,
    categoricalAsBooleansMap: null,
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
  switch (action.type) {
    /**********************************
        Keep a copy of 'universe'
  ***********************************/
    case "initialize success":
      return Object.assign({}, state, {
        _ranges: action.data.data.ranges,
        allGeneNames: action.data.data.genes
      });
    case "request cells success":
      const graphMap = {};
      const currentCellSelection = action.data.data.metadata.slice(0);
      _.each(action.data.data.graph, g => {
        graphMap[g[0]] = [g[1], g[2]];
      });

      /*
      construct a copy of the ranges object that only has categorical
      replace all counts with bool flags
      ie., everything starts out checked
      we mutate this map in the actions below
    */
      const categoricalAsBooleansMap = {},
        categoricalAsCellsMap = {};
      const continuousUserDefinedRanges = {};
      _.each(action.data.data.ranges, (value, key) => {
        if (
          key !== "CellName" &&
          value.options /* it's categorical, it has options instead of ranges */
        ) {
          const optionsAsBooleans = {},
            optionsAsCells = {};
          _.each(value.options, (_value, _key) => {
            optionsAsBooleans[_key] = true;
            optionsAsCells[_key] = [];
          });
          categoricalAsBooleansMap[key] = optionsAsBooleans;
          categoricalAsCellsMap[key] = optionsAsCells;
        } else if (key !== "CellName" && value.range) {
          continuousUserDefinedRanges[key] = null;
        }
      });

      _.each(currentCellSelection, cell => {
        cell["__selected__"] = true;
        cell["__color__"] =
          "rgba(0,0,0,1)"; /* initial color for all cells in all charts */
        cell["__colorRGB__"] = parseRGB(cell["__color__"]);

        // Add each cell to its categorical metadata set.
        _.forEach(cell, (_value, key) => {
          if (
            categoricalAsCellsMap[key] &&
            categoricalAsCellsMap[key][_value]
          ) {
            const s = categoricalAsCellsMap[key][_value];
            if (s) s.push(cell);
          }
        });
      });

      return Object.assign({}, state, {
        allCellsOnClient: action.data.data,
        currentCellSelection,
        graphMap,
        categoricalAsBooleansMap,
        categoricalAsCellsMap,
        continuousUserDefinedRanges,
        graphBrushSelection: null /* if we are getting new cells from the server, the layout (probably? definitely?) just changed, so this is now irrelevant, and we WILL need to call a function to reset state of this kind when cells success happens */
      });
    /* * * * * * * * * * * * * * * * * *
            User events
  * * * * * * * * * * * * * * * * * */
    case "parallel coordinates axes have been drawn":
      return Object.assign({}, state, {
        axesHaveBeenDrawn: true
      });
    case "continuous selection using parallel coords brushing": {
      return Object.assign({}, state, {
        continuousSelection: action.data,
        currentCellSelection:
          action.newSelection /* this comes from middleware */
      });
    }
    case "graph brush selection change":
      return Object.assign({}, state, {
        graphBrushSelection:
          action.brushCoords /* this has already been applied in middleware but store it for next time */,
        currentCellSelection:
          action.newSelection /* this comes from middleware */
      });
    case "graph brush deselect":
      return Object.assign({}, state, {
        graphBrushSelection: null,
        currentCellSelection:
          action.newSelection /* this comes from middleware */
      });
    case "continuous metadata histogram brush":
      return Object.assign({}, state, {
        newContinuousUserDefinedRanges:
          action.newContinuousUserDefinedRanges /* this has already been applied in middleware but store it for next time */,
        currentCellSelection:
          action.newSelection /* this comes from middleware */
      });
    case "change opacity deselected cells in 2d graph background":
      return Object.assign({}, state, {
        opacityForDeselectedCells: action.data
      });
    /*******************************
        Categorical metadata
  *******************************/
    case "categorical metadata filter select":
      return Object.assign({}, state, {
        categoricalAsBooleansMap:
          action.newCategoricalAsBooleansMap /* this comes from middleware */,
        currentCellSelection:
          action.newSelection /* this comes from middleware */
      });
    case "categorical metadata filter deselect":
      return Object.assign({}, state, {
        categoricalAsBooleansMap:
          action.newCategoricalAsBooleansMap /* this comes from middleware */,
        currentCellSelection:
          action.newSelection /* this comes from middleware */
      });
    case "categorical metadata filter none of these":
      return Object.assign({}, state, {
        categoricalAsBooleansMap:
          action.newCategoricalAsBooleansMap /* this comes from middleware */,
        currentCellSelection:
          action.newSelection /* this comes from middleware */
      });
    case "categorical metadata filter all of these":
      return Object.assign({}, state, {
        categoricalAsBooleansMap:
          action.newCategoricalAsBooleansMap /* this comes from middleware */,
        currentCellSelection:
          action.newSelection /* this comes from middleware */
      });
    /*******************************
            Color Scale
  *******************************/
    case "color by continuous metadata":
      return Object.assign({}, state, {
        colorAccessor: action.colorAccessor,
        currentCellSelection:
          action.currentSelectionWithUpdatedColors /* this comes from middleware */,
        colorScale: action.colorScale
      });
    case "color by expression":
      return Object.assign({}, state, {
        colorAccessor: action.gene,
        currentCellSelection:
          action.currentSelectionWithUpdatedColors /* this comes from middleware */,
        colorScale: action.colorScale
      });
    case "color by categorical metadata":
      return Object.assign({}, state, {
        colorAccessor:
          action.colorAccessor /* pass the scale through additionally, and it's a legend! */,
        currentCellSelection:
          action.currentSelectionWithUpdatedColors /* this comes from middleware */,
        colorScale: action.colorScale
      });
    case "store current cell selection as differential set 1":
      return Object.assign({}, state, {
        __storedStateForCelllist1__: action.data
      });
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
