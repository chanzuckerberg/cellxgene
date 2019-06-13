import * as d3 from "d3";
import styles from "./graph.css";

export default (responsive, graphPaddingRight, xy, text, colorBy) => {
  const containerWidth = responsive.width - graphPaddingRight;

  const svg = d3
    .select("#graphAttachPoint")
    .append("svg")
    .attr("id", "centroid-container")
    .attr("data-testid", "centroid-overlay")
    .attr("width", containerWidth)
    .attr("height", responsive.height)
    .attr("class", `${styles.graphSVG}`)
    .style("z-index", 998);
  //  TODO: Create own styles, ask Colin for an explanation on the css
  // For now I'm going to put centroid z-index at 998 and lasso on 999

  const label = svg
    .append("g")
    .attr("transform", `translate(${xy[0]}, ${xy[1]})`);

  label
    .append("text")
    .attr("text-anchor", "middle")
    .text(text)
    .style("font-family", "Roboto Condensed")
    .style("font-size", "18px")
    .style("font-weight", "700")
    .style("fill", colorBy ? "black" : "rgb(32, 178, 212)");

  return svg;
};
