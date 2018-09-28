// jshint esversion: 6
import * as d3 from "d3";
import styles from "./parallelCoordinates.css";
import { yAxis, brushstart } from "./util";

const drawAxes = (
  svg,
  ctx,
  dimensions,
  xscale,
  height,
  width,
  handleBrushAction,
  handleColorAction
) => {
  /*****************************************
  ******************************************
  Handles a brush event, toggling the display of foreground lines.
  ******************************************
  ******************************************/

  function brush() {
    const actives = [];
    svg
      .selectAll(".parcoords_axis .parcoords_brush")
      .filter(() => d3.brushSelection(this))
      .each(d => {
        actives.push({
          dimension: d,
          extent: d3.brushSelection(this)
        });
      });
    /* fire action, with selected dimensions & their values */
    handleBrushAction(actives);
  }

  const axes = svg
    .selectAll(".parcoords_axis")
    .data(dimensions)
    .enter()
    .append("g")
    .attr("class", `${styles.axis} parcoords_axis`)
    .attr("transform", (d, i) => `translate(${xscale(i)})`);

  axes
    .append("g")
    .each(d => {
      const renderAxis =
        "axis" in d
          ? d.axis.scale(d.scale) // custom axis
          : yAxis.scale(d.scale); // default axis
      d3.select(this).call(renderAxis);
    })
    .append("text")
    .on("click", d => {
      handleColorAction(d.key);
    })
    .attr("class", styles.title)
    .attr("text-anchor", "start")
    .text(d => ("description" in d ? `${d.description}  🖌️` : `${d.key}  🖌️`));

  // Add and store a brush for each axis.
  axes
    .append("g")
    .attr("class", `${styles.brush} parcoords_brush`)
    .each(d => {
      d3.select(this).call(
        (d.brush = d3
          .brushY()
          .extent([[-10, 0], [10, height]])
          .on("start", brushstart)
          .on("brush", brush)
          .on("end", brush))
      );
    })
    .selectAll("rect")
    .attr("x", -8)
    .attr("width", 16);

  return axes;
};

export default drawAxes;
