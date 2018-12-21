// jshint esversion: 6
import _ from "lodash";
import * as d3 from "d3";
import { interpolateRainbow, interpolateCool } from "d3-scale-chromatic";
import * as globals from "../globals";
import parseRGB from "../util/parseRGB";
import finiteExtent from "../util/finiteExtent";

/*
  https://medium.com/@jacobp100/you-arent-using-redux-middleware-enough-94ffe991e6
  storeInstance =>
    functionToCallWithAnActionThatWillSendItToTheNextMiddleware =>
    actionThatDispatchWasCalledWith =>
    valueToUseAsTheReturnValueOfTheDispatchCall
*/

/*
  What this file does:

  1. fire a filter action anywhere in the app
  2. ** this middleware checks to see the state of all the currently selected filters,
     including the new one
  3. ** create updated selection from a copy of all the cells presently on the client
     (this may be a subset of 'all')
  4. ** append that new selection to the action so that it magically appears in the reducer
     just because the action was fired

  This is nice because we keep a lot of filtering business logic centralized
  (what it means in practice to be selected)
*/

const updateCellColorsMiddleware = store => next => action => {
  const s = store.getState();

  /*
  this is a hardcoded map of the things we need to keep an eye on and update
  global cell selection in response to
  */
  const filterJustChanged =
    action.type === "color by expression" ||
    action.type === "color by continuous metadata" ||
    action.type === "color by categorical metadata";

  if (!filterJustChanged || !s.controls.world.obsAnnotations) {
    return next(
      action
    ); /* if the cells haven't loaded or the action wasn't a color change, bail */
  }

  const { obsAnnotations } = s.controls.world;
  let colorScale;
  const colorsByName = new Array(obsAnnotations.length);
  const colorsByRGB = new Array(obsAnnotations.length);

  /*
  in plain language...
   (a) once the cells have loaded.
   (b) each time a user changes a color control we need to update cellsMetadata colors
  This is available to all the draw functions as world.colorName[index] or world.colorRGB[index]
  */

  if (action.type === "color by categorical metadata") {
    const { categories } = _.filter(s.controls.world.schema.annotations.obs, {
      name: action.colorAccessor
    })[0];

    colorScale = d3
      .scaleSequential(interpolateRainbow)
      .domain([0, categories.length]);

    for (let i = 0; i < obsAnnotations.length; i += 1) {
      const obs = obsAnnotations[i];
      const c = colorScale(categories.indexOf(obs[action.colorAccessor]));
      colorsByName[i] = c;
      colorsByRGB[i] = parseRGB(c);
    }
  }

  if (action.type === "color by continuous metadata") {
    colorScale = d3
      .scaleLinear()
      .domain([0, action.rangeMaxForColorAccessor])
      .range([1, 0]);

    for (let i = 0; i < obsAnnotations.length; i += 1) {
      const val = obsAnnotations[i][action.colorAccessor];
      if (Number.isFinite(val)) {
        colorsByName[i] = interpolateCool(colorScale(val));
      } else {
        colorsByName[i] = globals.nonFiniteCellColor;
      }
      colorsByRGB[i] = parseRGB(colorsByName[i]);
    }
  }

  if (action.type === "color by expression") {
    const { gene, data } = action;
    const expression = data[gene]; // Float32Array
    colorScale = d3
      .scaleLinear()
      .domain(finiteExtent(expression))
      .range([
        1,
        0
      ]); /* invert viridis... probably pass this scale through to others */

    for (let i = 0, len = expression.length; i < len; i += 1) {
      const e = expression[i];
      if (Number.isFinite(e)) {
        colorsByName[i] = interpolateCool(colorScale(e));
      } else {
        colorsByName[i] = globals.nonFiniteCellColor;
      }
      colorsByRGB[i] = parseRGB(colorsByName[i]);
    }
  }

  /*
  append the result of all the filters to the action the user just triggered
  */
  const modifiedAction = Object.assign({}, action, {
    colors: { name: colorsByName, rgb: colorsByRGB },
    colorScale
  });

  return next(modifiedAction);
};

export default updateCellColorsMiddleware;
