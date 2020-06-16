/*
Helper functions for the embedded graph colors
*/
import * as d3 from "d3";
import { interpolateRainbow, interpolateCool } from "d3-scale-chromatic";
import memoize from "memoize-one";
import * as globals from "../../globals";
import parseRGB from "../parseRGB";
import finiteExtent from "../finiteExtent";
import { range } from "../range";

/*
given a color mode & accessor, generate an annoMatrix query that will
fulfill it
*/
export function createColorQuery(colorMode, colorByAccessor, schema) {
  if (!colorMode || !colorByAccessor || !schema) return null;
  switch (colorMode) {
    case "color by categorical metadata":
    case "color by continuous metadata": {
      return ["obs", colorByAccessor];
    }
    case "color by expression": {
      const varIndex = schema?.annotations?.var?.index;
      if (!varIndex) return null;
      return [
        "X",
        {
          field: "var",
          column: varIndex,
          value: colorByAccessor,
        },
      ];
    }
    default: {
      return null;
    }
  }
}

/*
create colors scale and RGB array and return as object. Parameters:
  * colorMode - categorical, etc.
  * colorByAccessor - the annotation label name
  * colorByDataframe - the actual color-by data
  * schema - the entire schema
  * userColors - optional user color table
Returns:
  {
    scale: color scale
    rgb: cell to color mapping
  }
*/
function _createColorTable(
  colorMode,
  colorByAccessor,
  colorByData,
  schema,
  userColors = null
) {
  switch (colorMode) {
    case "color by categorical metadata": {
      if (userColors && colorByAccessor in userColors) {
        return createUserColors(colorByData, colorByAccessor, userColors);
      }
      return createColorsByCategoricalMetadata(
        colorByData,
        colorByAccessor,
        schema
      );
    }
    case "color by continuous metadata": {
      return createColorsByContinuousMetadata(colorByData, colorByAccessor);
    }
    case "color by expression": {
      return createColorsByExpression(colorByData);
    }
    default: {
      const defaultCellColor = parseRGB(globals.defaultCellColor);
      return {
        rgb: new Array(schema.dataframe.nObs).fill(defaultCellColor),
        scale: undefined,
      };
    }
  }
}

export const createColorTable = memoize(_createColorTable);

export function loadUserColorConfig(userColors) {
  const convertedUserColors = {};
  Object.keys(userColors).forEach((category) => {
    const [colors, scaleMap] = Object.keys(userColors[category]).reduce(
      (acc, label, i) => {
        const color = parseRGB(userColors[category][label]);
        acc[0][label] = color;
        acc[1][i] = d3.rgb(255 * color[0], 255 * color[1], 255 * color[2]);
        return acc;
      },
      [{}, {}]
    );
    const scale = (i) => scaleMap[i];
    convertedUserColors[category] = { colors, scale };
  });
  return convertedUserColors;
}

function createUserColors(dataframe, colorAccessor, userColors) {
  const { colors, scale } = userColors[colorAccessor];
  const rgb = createRgbArray(dataframe, colors, colorAccessor);
  return { rgb, scale };
}

function createColorsByCategoricalMetadata(dataframe, colorAccessor, schema) {
  const { categories } = schema.annotations.obsByName[colorAccessor];

  const scale = d3
    .scaleSequential(interpolateRainbow)
    .domain([0, categories.length]);

  /* pre-create colors - much faster than doing it for each obs */
  const colors = categories.reduce((acc, cat, idx) => {
    acc[cat] = parseRGB(scale(idx));
    return acc;
  }, {});

  const rgb = createRgbArray(dataframe, colors, colorAccessor);
  return { rgb, scale };
}

export function createRgbArray(dataframe, colors, colorAccessor) {
  const data = dataframe.col(colorAccessor).asArray();
  const rgb = new Array(data.length);
  for (let i = 0, len = data.length; i < len; i += 1) {
    const label = data[i];
    rgb[i] = colors[label];
  }
  return rgb;
}

function createColorsByContinuousMetadata(dataframe, accessor) {
  const colorBins = 100;
  const col = dataframe.col(accessor);
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
  const data = col.asArray();
  const rgb = new Array(data.length);
  for (let i = 0, len = data.length; i < len; i += 1) {
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

function createColorsByExpression(dataframe) {
  const expression = dataframe.icol(0).asArray();
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

  const rgb = new Array(expression.length);
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
