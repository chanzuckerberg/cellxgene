import * as d3 from "d3";
import styles from "./graph.css";

export default (responsive, graphPaddingRight, xy, text) => {
  const containerWidth = responsive.width - graphPaddingRight;

  const svg = d3
    .select("#graphAttachPoint")
    .append("svg")
    .attr("id", "centroid")
    .attr("data-testid", "centroid-overlay")
    .attr("width", containerWidth)
    .attr("height", responsive.height)
    .attr("class", `${styles.graphSVG}`)
    .style("z-index", 998);
  //  TODO: Create own styles, ask Colin for an explanation on the css
  // For now I'm going to put centroid z-index at 998 and lasso on 999

  const radius = xy[2] * containerWidth * 1.5;

  const label = svg
    .append("g")
    .attr("transform", `translate(${xy[0]}, ${xy[1]})`);

  label
    .append("circle")
    .attr("r", radius)
    .style("fill", "rgba(244, 66, 66, .55)")
    .style("stroke", "rgba(0, 0, 0, .4)");

  label
    .append("text")
    .attr("text-anchor", "middle")
    .text(text)
    .style("font-size", "25px")
    .style("fill", "white")
    .style("stroke", "black");

  return svg;
};
