import _ from "lodash";
import * as d3 from "d3";

const paddingRight = 120;
const continuousChartWidth = 1200;

export const margin = { top: 66, right: 110, bottom: 20, left: 60 };
export const width =
  continuousChartWidth - margin.left - margin.right - paddingRight;
export const height = 340 - margin.top - margin.bottom;
export const innerHeight = height - 2;

export const devicePixelRatio = window.devicePixelRatio || 1;

export const createDimensions = (data) => {
  const newArr = [];
  _.each(data, (value, key) => {
    if (value.range) {
      newArr.push({
        key /* room for confusion: lodash calls this key, it's also the name of the property parallel coords code is looking for */,
        type: {
          within: (d, extent, dim) =>
            extent[0] <= dim.scale(d) && dim.scale(d) <= extent[1],
        },
        scale: d3
          .scaleSqrt()
          .range([innerHeight, 0])
          .domain([0, value.range.max]),
      });
    }
  });
  return newArr;
};

export const yAxis = d3.axisLeft();

export const brushstart = () => {
  d3.event.sourceEvent.stopPropagation();
};

// Unused.
// export const d3_functor = v => (typeof v === "function" ? v : () => v);
export const project = (d, dimensions, xscale) =>
  dimensions.map((p, i) => {
    // check if data element has property and contains a value
    if (!(p.key in d) || d[p.key] === null) return null;

    return [xscale(i), p.scale(d[p.key])];
  });
