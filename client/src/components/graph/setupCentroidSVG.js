import * as d3 from "d3";
import styles from "./graph.css";

export default (responsive, graphPaddingRight, labels, colorBy) => {
  const containerWidth = responsive.width - graphPaddingRight;

  const svg = d3
    .select("#graphAttachPoint")
    .append("svg")
    .attr("id", "centroid-container")
    .attr("data-testid", "centroid-overlay")
    .attr("width", containerWidth)
    .attr("height", responsive.height)
    .attr("class", `${styles.graphSVG}`)
    .style("z-index", 998)
    .style("pointer-events", "none");
  //  TODO: Create own styles, ask Colin for an explanation on the css
  // For now I'm going to put centroid z-index at 998 and lasso on 999

  // for (let i = 0, { length } = labels; i < length; i += 1) {
  const iter = labels.entries();
  let pair = iter.next().value;
  let value;
  let key;
  while (pair) {
    key = pair[0];
    value = pair[1];
    const label = svg
      .append("g")
      .attr("transform", `translate(${value[3]}, ${value[4]})`);

    label
      .append("text")
      .attr("text-anchor", "middle")
      .text(key)
      .style("font-family", "Roboto Condensed")
      .style("font-size", "18px")
      .style("font-weight", "700")
      .style("fill", colorBy ? "black" : "rgb(32, 178, 212)");
    pair = iter.next().value;
  }

  return svg;
};
