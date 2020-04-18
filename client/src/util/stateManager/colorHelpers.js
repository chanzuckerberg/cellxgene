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
  - colorMode - color-by mode. One of {null, "color by expression", "color by continuous metadata",
    "color by categorical metadata"}
  - colorAccessor - the obs annotations used for color-by
*/
export function createColors(world, colorMode = null, colorAccessor = null, userColors = null) {
  switch (colorMode) {
    case "color by categorical metadata": {
      if (userColors && colorAccessor in userColors) {
        return createUserColors(world, colorAccessor, userColors);
      }
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

export function loadUserColorConfig(userColors) {
  const convertedUserColors = {};
  Object.keys(userColors).forEach(category => {
    const [colors, scaleMap] = Object.keys(userColors[category]).reduce((acc, label, i) => {
      const color = parseRGB(userColors[category][label]);
      acc[0][label] = color;
      acc[1][i] = d3.rgb(255 * color[0], 255 * color[1], 255 * color[2]);
      return acc;
    }, [{}, {}]);
    const scale = i => scaleMap[i];
    convertedUserColors[category] = { colors, scale };
  });
  return convertedUserColors;
}

function createUserColors(world, colorAccessor, userColors) {
  const { colors, scale } = userColors[colorAccessor];
  const rgb = createRgbArray(world, colors, colorAccessor);
  return { rgb, scale };
}

function createColorsByCategoricalMetadata(world, colorAccessor) {
  const { categories } = world.schema.annotations.obsByName[colorAccessor];

  const scale = d3
    .scaleSequential(interpolateRainbow)
    .domain([0, categories.length]);

  /* pre-create colors - much faster than doing it for each obs */
  const colors = categories.reduce((acc, cat, idx) => {
    acc[cat] = parseRGB(scale(idx));
    return acc;
  }, {});

  const rgb = createRgbArray(world, colors, colorAccessor);
  return { rgb, scale };
}

export function createRgbArray(world, colors, colorAccessor) {
  const rgb = new Array(world.nObs);
  const df = world.obsAnnotations;
  const data = df.col(colorAccessor).asArray();
  for (let i = 0, len = df.length; i < len; i += 1) {
    const label = data[i];
    rgb[i] = colors[label];
  }
  return rgb;
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
