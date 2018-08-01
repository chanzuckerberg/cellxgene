// jshint esversion: 6
import styles from "./graph.css";
import _ from "lodash";
import * as globals from "../../globals";

/******************************************
*******************************************
          put svg & brush in DOM
*******************************************
******************************************/

export const setupSVGandBrushElements = (
  handleBrushSelectAction,
  handleBrushDeselectAction,
  responsive,
  graphPaddingTop
) => {
  const side = responsive.height - graphPaddingTop;
  const svg = d3
    .select("#graphAttachPoint")
    .append("svg")
    .attr("width", side)
    .attr("height", side)
    .attr("class", `${styles.graphSVG}`);

  const brush = d3
    .brush()
    .extent([
      [0, 0],
      [responsive.height - graphPaddingTop, responsive.height - graphPaddingTop]
    ])
    .on("brush", handleBrushSelectAction)
    .on("end", handleBrushDeselectAction);

  const brushContainer = svg
    .append("g")
    .attr("class", "graph_brush")
    .call(brush);

  return {
    svg,
    brushContainer,
    brush
  };
};
