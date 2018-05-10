// jshint esversion: 6
import styles from "./graph.css";
import _ from "lodash";
import * as globals from "../../globals";

/******************************************
*******************************************
          put svg in DOM
*******************************************
******************************************/

export const setupGraphElements = (
  handleBrushSelectAction,
  handleBrushDeselectAction
) => {
  var svg = d3
    .select("#graphAttachPoint")
    .append("svg")
    .attr("width", globals.graphWidth)
    .attr("height", globals.graphHeight)
    .attr("class", `${styles.graphSVG}`);

  setupGraphBrush(svg, handleBrushSelectAction, handleBrushDeselectAction);

  return {
    svg
  };
};

const setupGraphBrush = (
  svg,
  handleBrushSelectAction,
  handleBrushDeselectAction
) => {
  svg.append("g").call(
    d3
      .brush()
      .extent([[0, 0], [globals.graphWidth, globals.graphHeight]])
      .on("brush", handleBrushSelectAction)
      .on("end", handleBrushDeselectAction)
  );
};
