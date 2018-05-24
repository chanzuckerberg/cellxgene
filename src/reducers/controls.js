// jshint esversion: 6
import _ from "lodash";
import { parseRGB } from "../util/parseRGB";
var crossfilter = require("../util/typedCrossfilter");

// In the case where the REST server does not implement data schema
// declaration, we attempt to deduce it by sniffing the data.
//
function createSchemaByDataSniffing(ranges) {
  let schema = {};
  _.forEach(ranges, (value, key) => {
    schema[key] = {
      displayname: key,
      variabletype: value.options ? "categorical" : "continuous"
    };

    // Metadata field type is inferred by sniffing the data. This has some risks.
    // Caveats:
    //  * Values have been converted to native JS objects by the JSON parser.
    //  * Lots of assumptions about he REST API behaving properly (eg, min/max
    //    are the same type, etc).
    let type;
    if (schema[key].variabletype === "continuous" && value.range) {
      // Use min/max as a proxy for all data.
      const min = value.range.min;
      const max = value.range.max;
      type =
        typeof min !== "number" || typeof max !== "number"
          ? "string"
          : Number.isSafeInteger(min) && Number.isSafeInteger(max)
            ? "int"
            : "float";
    } else {
      // use an option value as a proxy for all data
      const aVal = value.options[0];
      type =
        typeof aVal !== "number"
          ? "string"
          : Number.isSafeInteger(aVal) ? "int" : "float";
    }
    schema[key].type = type;
  });
  return schema;
}

// Deduce the correct crossfilter dimension type from a metadata
// schema description.
//
function deduceDimensionType(attributes, fieldName) {
  let dimensionType;
  if (attributes.type === "string") {
    dimensionType = "enum";
  } else if (attributes.type === "int") {
    dimensionType = Int32Array;
  } else if (attributes.type === "float") {
    dimensionType = Float32Array;
  } else {
    console.error(
      `Warning - REST API returned unknown metadata schema (${
        attributes.type
      }) for field ${fieldName}.`
    );
    // skip it - we don't know what to do with this type
  }
  return dimensionType;
}

const Controls = (
  state = {
    _ranges: null /* this comes from initialize, this is universe */,
    allGeneNames: null,
    allCellsOnClient: null /* this comes from cells endpoint, this is world */,
    allCellsMetadata: null /* this comes from user actions, all draw components use this, it is created by middleware */,
    crossfilter: null /* the current user selection state */,
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
    case "initialize success": {
      if (!action.data.data.schema) {
        console.error("Warning - REST API omitted schema description.");
      }
      return Object.assign({}, state, {
        _ranges: action.data.data.ranges,
        allGeneNames: action.data.data.genes,
        schema: action.data.data.schema
      });
    }
    case "request cells success": {
      const allCellsMetadata = action.data.data.metadata.slice(0);
      const allCellsMetadataMap = _.keyBy(allCellsMetadata, "CellName");

      /*
      construct a copy of the ranges object that only has categorical
      replace all counts with bool flags
      ie., everything starts out checked
      we mutate this map in the actions below
      */
      const categoricalAsBooleansMap = {};
      _.each(action.data.data.ranges, (value, key) => {
        if (
          key !== "CellName" &&
          value.options /* it's categorical, it has options instead of ranges */
        ) {
          const optionsAsBooleans = {};
          _.each(value.options, (_value, _key) => {
            optionsAsBooleans[_key] = true;
          });
          categoricalAsBooleansMap[key] = optionsAsBooleans;
        }
      });

      const graph = action.data.data.graph;
      _.each(allCellsMetadata, (cell, idx) => {
        cell.__cellIndex__ = idx;
        cell.__color__ =
          "rgba(0,0,0,1)"; /* initial color for all cells in all charts */
        cell.__colorRGB__ = parseRGB(cell.__color__);
        cell.__x__ = graph[idx][1];
        cell.__y__ = graph[idx][2];
      });

      // Build the selection crossfilter.
      //
      let cellsCrossfilter = crossfilter(allCellsMetadata);
      let cellsDimensionsMap = {};
      cellsDimensionsMap.x = cellsCrossfilter.dimension(
        r => r.__x__,
        Float32Array
      );
      cellsDimensionsMap.y = cellsCrossfilter.dimension(
        r => r.__y__,
        Float32Array
      );

      // Now walk the schema and make an appropriate dimension for each
      // metadata field.   This is a simplistic mapping, and could be
      // optmized to use smaller scalars (to save memory) or larger
      // floating point where precision is needed.
      //
      // If we don't have a schema (bad server!), fake it by inferring
      // important fields from the ranges element.
      //
      if (!state.schema) {
        state.schema = createSchemaByDataSniffing(action.data.data.ranges);
      }
      _.forEach(state.schema, (attributes, key) => {
        if (key !== "CellName") {
          const dimensionType = deduceDimensionType(attributes, key);
          if (dimensionType) {
            cellsDimensionsMap[key] = cellsCrossfilter.dimension(
              r => r[key],
              dimensionType
            );
          }
        }
      });
      return Object.assign({}, state, {
        // this is only used as a flag that data has loaded.  Could be
        // removed (other variables would suffice for the same test).
        allCellsOnClient: action.data.data,
        allCellsMetadata,
        allCellsMetadataMap,
        categoricalAsBooleansMap,
        crossfilter: {
          cells: cellsCrossfilter,
          dimensionMap: cellsDimensionsMap
        },
        graphBrushSelection: null /* if we are getting new cells from the server, the layout (probably? definitely?) just changed, so this is now irrelevant, and we WILL need to call a function to reset state of this kind when cells success happens */
      });
    }
    /* * * * * * * * * * * * * * * * * *
            User events
  * * * * * * * * * * * * * * * * * */
    case "parallel coordinates axes have been drawn": {
      return Object.assign({}, state, {
        axesHaveBeenDrawn: true
      });
    }
    case "continuous selection using parallel coords brushing": {
      return Object.assign({}, state, {
        continuousSelection: action.data,
        crossfilter: {
          ...state.crossfilter
        }
      });
    }
    case "graph brush selection change": {
      state.crossfilter.dimensionMap.x.filterRange([
        action.brushCoords.northwest[0],
        action.brushCoords.southeast[0]
      ]);
      state.crossfilter.dimensionMap.y.filterRange([
        action.brushCoords.southeast[1],
        action.brushCoords.northwest[1]
      ]);
      return Object.assign({}, state, {
        graphBrushSelection: action.brushCoords,
        crossfilter: {
          ...state.crossfilter
        }
      });
    }
    case "graph brush deselect": {
      state.crossfilter.dimensionMap.x.filterAll();
      state.crossfilter.dimensionMap.y.filterAll();
      return Object.assign({}, state, {
        graphBrushSelection: null,
        crossfilter: {
          ...state.crossfilter
        }
      });
    }
    case "continuous metadata histogram brush": {
      // action.selection: metadata name being selected
      // action.range: filter range, or null if deselected
      if (!action.range) {
        state.crossfilter.dimensionMap[action.selection].filterAll();
      } else {
        state.crossfilter.dimensionMap[action.selection].filterRange(
          action.range
        );
      }
      return Object.assign({}, state, {
        crossfilter: {
          ...state.crossfilter
        }
      });
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
      state.crossfilter.dimensionMap[action.metadataField].filterEnum(
        _.filter(
          _.map(
            newCategoricalAsBooleansMap[action.metadataField],
            (val, key) => (val ? key : false)
          )
        )
      );
      return Object.assign({}, state, {
        categoricalAsBooleansMap: newCategoricalAsBooleansMap,
        crossfilter: {
          ...state.crossfilter
        }
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
      state.crossfilter.dimensionMap[action.metadataField].filterEnum(
        _.filter(
          _.map(
            newCategoricalAsBooleansMap[action.metadataField],
            (val, key) => (val ? key : false)
          )
        )
      );
      return Object.assign({}, state, {
        categoricalAsBooleansMap: newCategoricalAsBooleansMap,
        crossfilter: {
          ...state.crossfilter
        }
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
      state.crossfilter.dimensionMap[action.metadataField].filterNone();
      return Object.assign({}, state, {
        categoricalAsBooleansMap: newCategoricalAsBooleansMap,
        crossfilter: {
          ...state.crossfilter
        }
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
      state.crossfilter.dimensionMap[action.metadataField].filterAll();
      return Object.assign({}, state, {
        categoricalAsBooleansMap: newCategoricalAsBooleansMap,
        crossfilter: {
          ...state.crossfilter
        }
      });
    }
    /*******************************
            Color Scale
  *******************************/
    case "color by continuous metadata":
      return Object.assign({}, state, {
        colorAccessor: action.colorAccessor,
        allCellsMetadata:
          action.allCellsMetadataWithUpdatedColors /* this comes from middleware */,
        colorScale: action.colorScale
      });
    case "color by expression":
      return Object.assign({}, state, {
        colorAccessor: action.gene,
        allCellsMetadata:
          action.allCellsMetadataWithUpdatedColors /* this comes from middleware */,
        colorScale: action.colorScale
      });
    case "color by categorical metadata":
      return Object.assign({}, state, {
        colorAccessor:
          action.colorAccessor /* pass the scale through additionally, and it's a legend! */,
        allCellsMetadata:
          action.allCellsMetadataWithUpdatedColors /* this comes from middleware */,
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
