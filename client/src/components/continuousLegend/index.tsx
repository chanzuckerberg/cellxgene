import React from "react";
import { connect } from "react-redux";
import * as d3 from "d3";
import { interpolateCool } from "d3-scale-chromatic";

import {
  createColorTable,
  createColorQuery,
} from "../../util/stateManager/colorHelpers";

// create continuous color legend
const continuous = (selectorId: any, colorScale: any, colorAccessor: any) => {
  const legendHeight = 200;
  const legendWidth = 80;
  const margin = { top: 10, right: 60, bottom: 10, left: 2 };

  const canvas = d3
    .select(selectorId)
    .style("height", `${legendHeight}px`)
    .style("width", `${legendWidth}px`)
    .append("canvas")
    .attr("height", legendHeight - margin.top - margin.bottom)
    .attr("width", 1)
    .style("height", `${legendHeight - margin.top - margin.bottom}px`)
    .style("width", `${legendWidth - margin.left - margin.right}px`)
    .style("position", "absolute")
    .style("top", `${margin.top + 1}px`)
    .style("left", `${margin.left + 1}px`)
    .style(
      "transform",
      "scale(1,-1)"
    ) /* flip it! dark is high value light is low.
    we flip the color scale as well [1, 0] instead of [0, 1] */
    .node();

  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const ctx = canvas.getContext("2d");

  const legendScale = d3
    .scaleLinear()
    .range([1, legendHeight - margin.top - margin.bottom])
    .domain([
      colorScale.domain()[1],
      colorScale.domain()[0],
    ]); /* we flip this to make viridis colors dark if high in the color scale */

  // image data hackery based on http://bl.ocks.org/mbostock/048d21cf747371b11884f75ad896e5a5
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const image = ctx.createImageData(1, legendHeight);
  d3.range(legendHeight).forEach((i) => {
    const c = d3.rgb(colorScale(legendScale.invert(i)));
    image.data[4 * i] = c.r;
    image.data[4 * i + 1] = c.g;
    image.data[4 * i + 2] = c.b;
    image.data[4 * i + 3] = 255;
  });
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  ctx.putImageData(image, 0, 0);

  // A simpler way to do the above, but possibly slower. keep in mind the legend
  // width is stretched because the width attr of the canvas is 1
  // See http://stackoverflow.com/questions/4899799/whats-the-best-way-to-set-a-single-pixel-in-an-html5-canvas
  /*
  d3.range(legendheight).forEach(function(i) {
    ctx.fillStyle = colorscale(legendscale.invert(i));
    ctx.fillRect(0,i,1,1);
  });
  */

  const legendAxis = d3
    .axisRight(legendScale)
    .ticks(6)
    .tickFormat(
      d3.format(
        legendScale.domain().some((n) => Math.abs(n) >= 10000) ? ".0e" : ","
      )
    );

  const svg = d3
    .select(selectorId)
    .append("svg")
    .attr("height", `${legendHeight}px`)
    .attr("width", `${legendWidth}px`)
    .style("position", "absolute")
    .style("left", "0px")
    .style("top", "0px");

  svg
    .append("g")
    .attr("class", "axis")
    .attr(
      "transform",
      `translate(${legendWidth - margin.left - margin.right + 3},${margin.top})`
    )
    .call(legendAxis);

  // text label for the y axis
  svg
    .append("text")
    .attr("transform", "rotate(-90)")
    .attr("y", 2)
    .attr("x", 0 - legendHeight / 2)
    .attr("dy", "1em")
    .attr("data-testid", "continuous_legend_color_by_label")
    .attr("aria-label", colorAccessor)
    .style("text-anchor", "middle")
    .style("fill", "white")
    .text(colorAccessor);
};

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  annoMatrix: (state as any).annoMatrix,
  colors: (state as any).colors,
  genesets: (state as any).genesets.genesets,
}))
class ContinuousLegend extends React.Component {
  async componentDidUpdate(prevProps: any) {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Read... Remove this comment to see the full error message
    const { annoMatrix, colors, genesets } = this.props;
    if (!colors || !annoMatrix) return;
    if (colors !== prevProps?.colors || annoMatrix !== prevProps?.annoMatrix) {
      const { schema } = annoMatrix;
      const { colorMode, colorAccessor, userColors } = colors;
      const colorQuery = createColorQuery(
        colorMode,
        colorAccessor,
        schema,
        genesets
      );
      const colorDf = colorQuery ? await annoMatrix.fetch(...colorQuery) : null;
      const colorTable = createColorTable(
        colorMode,
        colorAccessor,
        colorDf,
        schema,
        userColors
      );
      const colorScale = colorTable.scale;
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'range' does not exist on type '((idx: an... Remove this comment to see the full error message
      const range = colorScale?.range;
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'domain' does not exist on type '((idx: a... Remove this comment to see the full error message
      const [domainMin, domainMax] = colorScale?.domain?.() ?? [0, 0];
      /* always remove it, if it's not continuous we don't put it back. */
      d3.select("#continuous_legend").selectAll("*").remove();
      if (colorAccessor && colorScale && range && domainMin < domainMax) {
        /* fragile! continuous range is 0 to 1, not [#fa4b2c, ...], make this a flag? */
        if (range()[0][0] !== "#") {
          continuous(
            "#continuous_legend",
            // @ts-expect-error ts-migrate(2339) FIXME: Property 'domain' does not exist on type '((idx: a... Remove this comment to see the full error message
            d3.scaleSequential(interpolateCool).domain(colorScale.domain()),
            colorAccessor
          );
        }
      }
    }
  }

  render() {
    return (
      <div
        id="continuous_legend"
        style={{
          position: "absolute",
          left: 8,
          top: 35,
          zIndex: 1,
          pointerEvents: "none",
        }}
      />
    );
  }
}

export default ContinuousLegend;
