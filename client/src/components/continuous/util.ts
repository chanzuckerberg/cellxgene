import each from "lodash.foreach";
import * as d3 from "d3";

const paddingRight = 120;
const continuousChartWidth = 1200;

export const margin = { top: 66, right: 110, bottom: 20, left: 60 };
export const width =
  continuousChartWidth - margin.left - margin.right - paddingRight;
export const height = 340 - margin.top - margin.bottom;
export const innerHeight = height - 2;

export const devicePixelRatio = window.devicePixelRatio || 1;

export const createDimensions = (data: any) => {
  const newArr: any = [];
  each(data, (value, key) => {
    if (value.range) {
      newArr.push({
        key /* room for confusion: lodash calls this key, it's also the name of the property parallel coords code is looking for */,
        type: {
          within: (d: any, extent: any, dim: any) =>
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

// @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 0.
export const yAxis = d3.axisLeft();

export const brushstart = () => {
  (d3 as any).event.sourceEvent.stopPropagation();
};

// Unused.
// export const d3_functor = v => (typeof v === "function" ? v : () => v);
export const project = (d: any, dimensions: any, xscale: any) =>
  dimensions.map((p: any, i: any) => {
    // check if data element has property and contains a value
    if (!(p.key in d) || d[p.key] === null) return null;

    return [xscale(i), p.scale(d[p.key])];
  });
