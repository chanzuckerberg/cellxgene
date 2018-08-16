// jshint esversion: 6
import uri from "urijs";
import * as globals from "../globals";
import _ from "lodash";
import { parseRGB } from "../util/parseRGB";
import * as d3 from "d3";
import { interpolateViridis } from "d3-scale-chromatic";

/*
  https://medium.com/@jacobp100/you-arent-using-redux-middleware-enough-94ffe991e6
  storeInstance => functionToCallWithAnActionThatWillSendItToTheNextMiddleware => actionThatDispatchWasCalledWith => valueToUseAsTheReturnValueOfTheDispatchCall
*/

/*
  What this file does:

  1. fire a filter action anywhere in the app
  2. ** this middleware checks to see the state of all the currently selected filters, including the new one
  3. ** create updated selection from a copy of all the cells presently on the client (this may be a subset of 'all')
  4. ** append that new selection to the action so that it magically appears in the reducer just because the action was fired

  This is nice because we keep a lot of filtering business logic centralized (what it means in practice to be selected)
*/

const updateCellColorsMiddleware = store => {
  return next => {
    return action => {
      const s = store.getState();

      /* this is a hardcoded map of the things we need to keep an eye on and update global cell selection in response to */
      const filterJustChanged =
        action.type === "color by expression" ||
        action.type === "color by continuous metadata" ||
        action.type === "color by categorical metadata";

      if (!filterJustChanged || !s.controls.cellsMetadata) {
        return next(
          action
        ); /* if the cells haven't loaded or the action wasn't a color change, bail */
      }

      let cellsMetadataWithUpdatedColors = s.controls.cellsMetadata.slice(0);
      let colorScale;

      /*
         in plain language...

         (a) once the cells have loaded.
         (b) each time a user changes a color control we need to update cellsMetadata colors

         This is available to all the draw functions as cell["__color__"] and cell["__colorRGB__"]
      */

      if (action.type === "color by categorical metadata") {
        colorScale = d3.scaleOrdinal().range(globals.ordinalColors);

        for (let i = 0; i < cellsMetadataWithUpdatedColors.length; i++) {
          const cell = cellsMetadataWithUpdatedColors[i];
          let c = colorScale(cell[action.colorAccessor]);
          cell.__color__ = c;
          cell.__colorRGB__ = parseRGB(c);
        }
      }

      if (action.type === "color by continuous metadata") {
        colorScale = d3
          .scaleLinear()
          .domain([0, action.rangeMaxForColorAccessor])
          .range([1, 0]);

        _.each(cellsMetadataWithUpdatedColors, (cell, i) => {
          let c = interpolateViridis(colorScale(cell[action.colorAccessor]));
          cellsMetadataWithUpdatedColors[i]["__color__"] = c;
          cellsMetadataWithUpdatedColors[i]["__colorRGB__"] = parseRGB(c);
        });
      }

      if (action.type === "color by expression") {
        const indexOfGene = 0; /* we only get one, this comes from server as needed now */

        const expressionMap = {};
        /*
          converts [{cellname: cell123, e}, {}]

          expressionMap = {
            cell123: [123, 2],
            cell789: [0, 8]
          }
        */
        _.each(action.data.data.cells, cell => {
          /* this action is coming directly from the server */
          expressionMap[cell.cellname] = cell.e;
        });

        const minExpressionCell = _.minBy(action.data.data.cells, cell => {
          return cell.e[indexOfGene];
        });

        const maxExpressionCell = _.maxBy(action.data.data.cells, cell => {
          return cell.e[indexOfGene];
        });

        // console.log('middle', action, expressionMap, minExpressionCell)

        colorScale = d3
          .scaleLinear()
          .domain([
            minExpressionCell.e[indexOfGene],
            maxExpressionCell.e[indexOfGene]
          ])
          .range([
            1,
            0
          ]); /* invert viridis... probably pass this scale through to others */

        _.each(cellsMetadataWithUpdatedColors, (cell, i) => {
          let c = interpolateViridis(
            colorScale(expressionMap[cell.CellName][indexOfGene])
          );
          cellsMetadataWithUpdatedColors[i]["__color__"] = c;
          cellsMetadataWithUpdatedColors[i]["__colorRGB__"] = parseRGB(c);
        });
      }

      /*
        append the result of all the filters to the action the user just triggered
      */
      let modifiedAction = Object.assign({}, action, {
        cellsMetadataWithUpdatedColors,
        colorScale
      });

      return next(modifiedAction);
    };
  };
};

export default updateCellColorsMiddleware;
