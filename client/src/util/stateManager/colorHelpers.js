/*
Helper functions for the embedded graph colors
*/
import * as d3 from "d3";
import { interpolateRainbow, interpolateCool } from "d3-scale-chromatic";
import * as globals from "../../globals";
import parseRGB from "../parseRGB";
import finiteExtent from "../finiteExtent";
import { range } from "../range";

/*
create new colors state object.   Paramters:
  - world - current world object
  - mode - color-by mode.  One of: null, "color by expression",
    "color by continuous metadata", "color by categorical metadata"
  -
*/
export function createColors(world, colorMode = null, colorAccessor = null) {
  switch (colorMode) {
    case "color by categorical metadata": {
      return createColorsByCategoricalMetadata(world, colorAccessor);
    }
    case "color by continuous metadata": {
      return createColorsByContinuousMetadata(world, colorAccessor);
    }
    case "color by expression": {
      return createColorsByExpression(world, colorAccessor);
    }
    default: {
      const defaultCellColor = parseRGB(globals.defaultCellColor);
      return {
        rgb: new Array(world.nObs).fill(defaultCellColor),
        scale: undefined
      };
    }
  }
}

function createColorsByCategoricalMetadata(world, accessor) {
  const { categories } = world.schema.annotations.obsByName[accessor];

  const scale = d3
    .scaleSequential(interpolateRainbow)
    .domain([0, categories.length]);

  /* pre-create colors - much faster than doing it for each obs */
  const colors = categories.reduce((acc, cat, idx) => {
    acc[cat] = parseRGB(scale(idx));
    return acc;
  }, {});

  const rgb = new Array(world.nObs);
  const df = world.obsAnnotations;
  const data = df.col(accessor).asArray();
  for (let i = 0, len = df.length; i < len; i += 1) {
    const cat = data[i];
    rgb[i] = colors[cat];
  }
  return { rgb, scale };
}

function createColorsByContinuousMetadata(world, accessor) {
  const colorBins = 100;
  const col = world.obsAnnotations.col(accessor);
  const { min, max } = col.summarize();
  const scale = d3
    .scaleQuantile()
    .domain([min, max])
    .range(range(colorBins - 1, -1, -1));

  /* pre-create colors - much faster than doing it for each obs */
  const colors = new Array(colorBins);
  for (let i = 0; i < colorBins; i += 1) {
    colors[i] = parseRGB(interpolateCool(i / colorBins));
  }

  const nonFiniteColor = parseRGB(globals.nonFiniteCellColor);
  const rgb = new Array(world.nObs);
  const data = col.asArray();
  for (let i = 0, len = world.obsAnnotations.length; i < len; i += 1) {
    const val = data[i];
    if (Number.isFinite(val)) {
      const c = scale(val);
      rgb[i] = colors[c];
    } else {
      rgb[i] = nonFiniteColor;
    }
  }
  return { rgb, scale };
}

function createColorsByExpression(world, accessor) {
  const expression = world.varData.col(accessor).asArray();
  const colorBins = 100;
  const [min, max] = finiteExtent(expression);
  const scale = d3
    .scaleQuantile()
    .domain([min, max])
    .range(range(colorBins - 1, -1, -1));

  /* pre-create colors - much faster than doing it for each obs */
  const colors = new Array(colorBins);
  for (let i = 0; i < colorBins; i += 1) {
    colors[i] = parseRGB(interpolateCool(i / colorBins));
  }
  const nonFiniteColor = parseRGB(globals.nonFiniteCellColor);

  const rgb = new Array(world.nObs);
  for (let i = 0, len = expression.length; i < len; i += 1) {
    const e = expression[i];
    if (Number.isFinite(e)) {
      const c = scale(e);
      rgb[i] = colors[c];
    } else {
      rgb[i] = nonFiniteColor;
    }
  }
  return { rgb, scale };
}

export const resetColors = world => {
  const { rgb, scale } = createColors(world);
  return {
    colorMode: null,
    colorAccessor: null,
    rgb,
    scale
  };
};

export const checkIfColorByDiffexpAndResetColors = (
  prevControls,
  state,
  prevWorld
) => {
  if (prevControls.diffexpGenes.includes(state.colorAccessor)) {
    return {
      ...state,
      ...resetColors(prevWorld)
    };
  }
  return null;
};
