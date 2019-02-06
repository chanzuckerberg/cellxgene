// jshint esversion: 6
import * as d3 from "d3";
import styles from "./graph.css";
import Lasso from "./setupLasso";

/******************************************
*******************************************
          put svg & brush in DOM
*******************************************
******************************************/

export default (
  handleBrushSelectAction,
  handleBrushDeselectAction,
  responsive,
  graphPaddingRight,
  handleLassoEnd,
  handleLassoStart
) => {
  const svg = d3
    .select("#graphAttachPoint")
    .append("svg")
    .attr("width", responsive.width - graphPaddingRight)
    .attr("height", responsive.height)
    .attr("class", `${styles.graphSVG}`);

  const brush = d3
    .brush()
    .extent([[0, 0], [responsive.width - graphPaddingRight, responsive.height]])
    .on("brush", handleBrushSelectAction)
    .on("end", handleBrushDeselectAction);

  const brushContainer = svg
    .append("g")
    .attr("class", "graph_brush")
    .call(brush);

  const lassoInstance = Lasso()
    .on("end", handleLassoEnd)
    .on("start", handleLassoStart);

  const lasso = svg.call(lassoInstance);

  return {
    svg,
    brushContainer,
    brush,
    lasso
  };
};
