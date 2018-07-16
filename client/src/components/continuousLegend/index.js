// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";

// create continuous color legend
// http://bl.ocks.org/syntagmatic/e8ccca52559796be775553b467593a9f
const continuous = (selector_id, colorscale) => {
  var legendheight = 200,
    legendwidth = 80,
    margin = { top: 10, right: 60, bottom: 10, left: 2 };

  var canvas = d3
    .select(selector_id)
    .style("height", legendheight + "px")
    .style("width", legendwidth + "px")
    // .style("position", "relative")
    .append("canvas")
    .attr("height", legendheight - margin.top - margin.bottom)
    .attr("width", 1)
    .style("height", legendheight - margin.top - margin.bottom + "px")
    .style("width", legendwidth - margin.left - margin.right + "px")
    // .style("border", "1px solid #000")
    .style("position", "absolute")
    .style("top", margin.top + 1 + "px")
    .style("left", margin.left + 1 + "px")
    .style(
      "transform",
      "scale(1,-1)"
    ) /* flip it! dark is high value light is low. we flip the color scale as well [1, 0] instead of [0, 1] */
    .node();

  var ctx = canvas.getContext("2d");

  var legendscale = d3
    .scaleLinear()
    .range([1, legendheight - margin.top - margin.bottom])
    .domain([
      colorscale.domain()[1],
      colorscale.domain()[0]
    ]); /* we flip this to make viridis colors dark if high in the color scale */

  // image data hackery based on http://bl.ocks.org/mbostock/048d21cf747371b11884f75ad896e5a5
  var image = ctx.createImageData(1, legendheight);
  d3.range(legendheight).forEach(function(i) {
    var c = d3.rgb(colorscale(legendscale.invert(i)));
    image.data[4 * i] = c.r;
    image.data[4 * i + 1] = c.g;
    image.data[4 * i + 2] = c.b;
    image.data[4 * i + 3] = 255;
  });
  ctx.putImageData(image, 0, 0);

  // A simpler way to do the above, but possibly slower. keep in mind the legend width is stretched because the width attr of the canvas is 1
  // See http://stackoverflow.com/questions/4899799/whats-the-best-way-to-set-a-single-pixel-in-an-html5-canvas
  /*
  d3.range(legendheight).forEach(function(i) {
    ctx.fillStyle = colorscale(legendscale.invert(i));
    ctx.fillRect(0,i,1,1);
  });
  */

  var legendaxis = d3
    .axisRight()
    .scale(legendscale)
    .tickSize(6)
    .ticks(8);

  var svg = d3
    .select(selector_id)
    .append("svg")
    .attr("height", legendheight + "px")
    .attr("width", legendwidth + "px")
    .style("position", "absolute")
    .style("left", "0px")
    .style("top", "0px");

  svg
    .append("g")
    .attr("class", "axis")
    .attr(
      "transform",
      "translate(" +
        (legendwidth - margin.left - margin.right + 3) +
        "," +
        margin.top +
        ")"
    )
    .call(legendaxis);
};

@connect(state => {
  return {
    colorAccessor: state.controls.colorAccessor,
    colorScale: state.controls.colorScale,
    responsive: state.responsive
  };
})
class ContinuousLegend extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }
  componentWillReceiveProps(nextProps) {
    if (
      nextProps.colorAccessor !== this.props.colorAccessor ||
      nextProps.responsive.height !== this.props.responsive.height ||
      nextProps.responsive.width !== this.props.responsive.width
    ) {
      /* always remove it, if it's not continuous we don't put it back. */
      d3
        .select("#continuous_legend")
        .selectAll("*")
        .remove();
    }

    if (nextProps.colorAccessor && nextProps.colorScale) {
      /* fragile! continuous range is 0 to 1, not [#fa4b2c, ...], make this a flag? */
      if (nextProps.colorScale.range()[0][0] !== "#") {
        continuous(
          "#continuous_legend",
          d3
            .scaleSequential(d3.interpolateViridis)
            .domain(nextProps.colorScale.domain())
        );
      }
    }
  }
  drawScale() {}
  render() {
    return (
      <div
        id="continuous_legend"
        style={{
          position: "fixed",
          display: this.props.colorAccessor ? "inherit" : "none",
          right: 0,
          top: this.props.responsive.height / 2
        }}
      />
    );
  }
}

export default ContinuousLegend;
