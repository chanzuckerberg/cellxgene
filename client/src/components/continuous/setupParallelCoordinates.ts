/*****************************************
******************************************
      Setup SVG & Canvas elements
******************************************
******************************************/
import * as d3 from "d3";

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
const setupParallelCoordinates = (width: any, height: any, margin: any) => {
  const container = d3.select("#parcoords");

  const svg = container
    .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);

  const canvas = container
    .append("canvas")
    .attr("width", width * devicePixelRatio)
    .attr("height", height * devicePixelRatio)
    .style("width", `${width}px`)
    .style("height", `${height}px`)
    .style("margin-top", `${margin.top}px`)
    .style("margin-left", `${margin.left}px`);

  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const ctx = canvas.node().getContext("2d");
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  ctx.globalCompositeOperation = "darken";
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  ctx.globalAlpha = 0.15;
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  ctx.lineWidth = 1.5;
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  ctx.scale(devicePixelRatio, devicePixelRatio);

  return {
    svg,
    ctx,
  };
};

export default setupParallelCoordinates;
