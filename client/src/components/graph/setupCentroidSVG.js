import * as d3 from "d3";
import styles from "./graph.css";

export default (responsive, graphPaddingRight, labels) => {
  const containerWidth = responsive.width - graphPaddingRight;

  const svg = d3
    .select("#graphAttachPoint")
    .select("svg")
    .attr("width", containerWidth)
    .attr("height", responsive.height)
    .attr("class", `${styles.graphSVG}`)
    .style("background", "rgba(255, 255, 255, 0.625");
  //  TODO: Create own styles, ask Colin for an explanation on the css
  // For now I'm going to put centroid z-index at 998 and lasso on 999

  // Iterate over the categoryValue -> coordinates, Map
  const iter = labels.entries();
  let pair = iter.next().value;
  let value;
  let key;
  while (pair) {
    key = pair[0];
    value = pair[1];
    // Create the label using the screen calculated coordinates
    const label = svg
      .append("g")
      .attr("class", "centroid-label")
      .attr("transform", `translate(${value[2]}, ${value[3]})`)
    // And the key as the text
    label
      .append("text")
      .attr("text-anchor", "middle")
      .attr("id", `svg${key.replace(/[^\w]/gi, "")}-label`)
      .text(key.length > 20 ? `${key.substr(0, 20)}...` : key)
      .style("font-family", "Roboto Condensed")
      .style("font-size", "12px")
      .style("fill", "black");
    pair = iter.next().value;
  }

  return svg;
};
