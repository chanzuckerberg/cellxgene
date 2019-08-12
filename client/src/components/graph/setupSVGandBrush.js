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
  selectionToolType,
  handleStartAction,
  handleDragAction,
  handleEndAction,
  handleCancelAction,
  responsive,
  graphPaddingRight,
  graphInteractionMode
) => {
  const svg = d3
    .select("#graphAttachPoint")
    .select("svg")
    .attr("data-testid", "layout-overlay")
    .attr("width", responsive.width - graphPaddingRight)
    .attr("height", responsive.height)
    .attr("class", `${styles.graphSVG}`)
    .style("z-index", 999)
    .attr(
      "pointer-events",
      graphInteractionMode === "select" ? "auto" : "none"
    );
  /* .style("display", graphInteractionMode === "select" ? "inherit" : "none"); */

  if (selectionToolType === "brush") {
    const brush = d3
      .brush()
      .extent([
        [0, 0],
        [responsive.width - graphPaddingRight, responsive.height]
      ])
      .on("start", handleStartAction)
      .on("brush", handleDragAction)
      // FYI, brush doesn't generate cancel
      .on("end", handleEndAction);

    const brushContainer = svg
      .append("g")
      .attr("class", "graph_brush")
      .call(brush);

    return { svg, container: brushContainer, tool: brush };
  }

  if (selectionToolType === "lasso") {
    const lasso = Lasso()
      .on("end", handleEndAction)
      // FYI, Lasso doesn't generate drag
      .on("start", handleStartAction)
      .on("cancel", handleCancelAction);

    const lassoContainer = svg.call(lasso);

    return { svg, container: lassoContainer, tool: lasso };
  }

  throw new Error("unknown graph selection tool");
};
