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

  const obsAnnotations = _.get(s.controls, "world.obsAnnotations", null);
  if (!filterJustChanged || !obsAnnotations) {
    return next(
      action
    ); /* if the cells haven't loaded or the action wasn't a color change, bail */
  }

  let colorScale;
  const colorsByRGB = new Array(obsAnnotations.length);

  /*
  in plain language...
   (a) once the cells have loaded.
   (b) each time a user changes a color control we need to update cellsMetadata colors
  This is available to all the draw functions as controls.colorRGB[index]
  */

  if (action.type === "color by categorical metadata") {
    const { categories } = _.filter(s.controls.world.schema.annotations.obs, {
      name: action.colorAccessor
    })[0];

    colorScale = d3
      .scaleSequential(interpolateRainbow)
      .domain([0, categories.length]);

    /* pre-create colors - much faster than doing it for each obs */
    const colors = _.transform(categories, (acc, cat, idx) => {
      acc[cat] = parseRGB(colorScale(idx));
    });

    const key = action.colorAccessor;
    const col = obsAnnotations.col(key).asArray();
    for (let i = 0, len = obsAnnotations.length; i < len; i += 1) {
      const cat = col[i];
      colorsByRGB[i] = colors[cat];
    }
  }

  if (action.type === "color by continuous metadata") {
    const colorBins = 100;
    const [min, max] = [0, action.rangeMaxForColorAccessor];
    colorScale = d3
      .scaleQuantile()
      .domain([min, max])
      .range(_.range(colorBins - 1, -1, -1));

    /* pre-create colors - much faster than doing it for each obs */
    const colors = new Array(colorBins);
    for (let i = 0; i < colorBins; i += 1) {
      colors[i] = parseRGB(interpolateCool(i / colorBins));
    }

    const key = action.colorAccessor;
    const nonFiniteColor = parseRGB(globals.nonFiniteCellColor);
    const col = obsAnnotations.col(key).asArray();
    for (let i = 0, len = obsAnnotations.length; i < len; i += 1) {
      const val = col[i];
      if (Number.isFinite(val)) {
        const c = colorScale(val);
        colorsByRGB[i] = colors[c];
      } else {
        colorsByRGB[i] = nonFiniteColor;
      }
    }
  }

  if (action.type === "color by expression") {
    const { gene, data } = action;
    const expression = data[gene]; // Float32Array
    const colorBins = 100;
    const [min, max] = finiteExtent(expression);
    colorScale = d3
      .scaleQuantile()
      .domain([min, max])
      .range(_.range(colorBins - 1, -1, -1));

    /* pre-create colors - much faster than doing it for each obs */
    const colors = new Array(colorBins);
    for (let i = 0; i < colorBins; i += 1) {
      colors[i] = parseRGB(interpolateCool(i / colorBins));
    }
    const nonFiniteColor = parseRGB(globals.nonFiniteCellColor);

    for (let i = 0, len = expression.length; i < len; i += 1) {
      const e = expression[i];
      if (Number.isFinite(e)) {
        const c = colorScale(e);
        colorsByRGB[i] = colors[c];
      } else {
        colorsByRGB[i] = nonFiniteColor;
      }
    }
  }

  /*
  append the result of all the filters to the action the user just triggered
  */
  const modifiedAction = Object.assign({}, action, {
    colors: { rgb: colorsByRGB },
    colorScale
  });

  return next(modifiedAction);
};

export default updateCellColorsMiddleware;
