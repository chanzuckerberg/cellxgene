// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import * as d3 from "d3";
import { interpolateViridis } from "d3-scale-chromatic";

// create continuous color legend
// http://bl.ocks.org/syntagmatic/e8ccca52559796be775553b467593a9f
const continuous = (selectorId, colorscale) => {
  const legendheight = 200;
  const legendwidth = 80;
  const margin = { top: 10, right: 60, bottom: 10, left: 2 };

  const canvas = d3
    .select(selectorId)
    .style("height", `${legendheight}px`)
    .style("width", `${legendwidth}px`)
    .append("canvas")
    .attr("height", legendheight - margin.top - margin.bottom)
    .attr("width", 1)
    .style("height", `${legendheight - margin.top - margin.bottom}px`)
    .style("width", `${legendwidth - margin.left - margin.right}px`)
    .style("position", "absolute")
    .style("top", `${margin.top + 1}px`)
    .style("left", `${margin.left + 1}px`)
    .style(
      "transform",
      "scale(1,-1)"
    ) /* flip it! dark is high value light is low.
    we flip the color scale as well [1, 0] instead of [0, 1] */
    .node();

  const ctx = canvas.getContext("2d");

  const legendscale = d3
    .scaleLinear()
    .range([1, legendheight - margin.top - margin.bottom])
    .domain([
      colorscale.domain()[1],
      colorscale.domain()[0]
    ]); /* we flip this to make viridis colors dark if high in the color scale */

  // image data hackery based on http://bl.ocks.org/mbostock/048d21cf747371b11884f75ad896e5a5
  const image = ctx.createImageData(1, legendheight);
  d3.range(legendheight).forEach(i => {
    const c = d3.rgb(colorscale(legendscale.invert(i)));
    image.data[4 * i] = c.r;
    image.data[4 * i + 1] = c.g;
    image.data[4 * i + 2] = c.b;
    image.data[4 * i + 3] = 255;
  });
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

  const legendaxis = d3
    .axisRight()
    .scale(legendscale)
    .tickSize(6)
    .ticks(8);

  const svg = d3
    .select(selectorId)
    .append("svg")
    .attr("height", `${legendheight}px`)
    .attr("width", `${legendwidth}px`)
    .style("position", "absolute")
    .style("left", "0px")
    .style("top", "0px");

  svg
    .append("g")
    .attr("class", "axis")
    .attr(
      "transform",
      `translate(${legendwidth - margin.left - margin.right + 3},${margin.top})`
    )
    .call(legendaxis);
};

@connect(state => ({
  colorAccessor: state.controls.colorAccessor,
  colorScale: state.controls.colorScale,
  responsive: state.responsive
}))
class ContinuousLegend extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  componentDidUpdate(prevProps) {
    const { colorAccessor, responsive, colorScale } = this.props;
    if (
      prevProps.colorAccessor !== colorAccessor ||
      prevProps.responsive.height !== responsive.height ||
      prevProps.responsive.width !== responsive.width
    ) {
      /* always remove it, if it's not continuous we don't put it back. */
      d3.select("#continuous_legend")
        .selectAll("*")
        .remove();
    }

    if (colorAccessor && colorScale) {
      /* fragile! continuous range is 0 to 1, not [#fa4b2c, ...], make this a flag? */
      if (colorScale.range()[0][0] !== "#") {
        continuous(
          "#continuous_legend",
          d3.scaleSequential(interpolateViridis).domain(colorScale.domain())
        );
      }
    }
  }

  render() {
    const { colorAccessor, responsive } = this.props;
    return (
      <div
        id="continuous_legend"
        style={{
          position: "fixed",
          display: colorAccessor ? "inherit" : "none",
          right: 0,
          top: responsive.height / 2
        }}
      />
    );
  }
}

export default ContinuousLegend;
