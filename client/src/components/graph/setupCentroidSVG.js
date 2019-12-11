import * as d3 from "d3";
import styles from "./graph.css";

export default (
  responsive,
  graphPaddingRight,
  labels,
  handleMouseEnter,
  handleMouseExit
) => {
  const centroids = [];
  const fontSize = "20px";
  const container = "#model-transformation-group";
  // Iterate over the categoryValue -> coordinates, Map
  const iter = labels.entries();
  let pair = iter.next().value;
  let value;
  let key;
  while (pair) {
    key = pair[0];
    value = pair[1];
    // Create the label using the screen calculated coordinates
    centroids.push(
      d3
        .select(container)
        .append("g")
        .attr("class", "centroid-label")
        .attr("transform", `translate(${value[0]}, ${value[1]})`)
        .datum({ key })
        .on("mouseover", handleMouseEnter)
        .on("mouseout", handleMouseExit)
        .append("text") // And the key as the text
        .attr("text-anchor", "middle")
        .attr("id", `svg${key.replace(/[^\w]/gi, "")}-label`)
        .text(key.length > 20 ? `${key.substr(0, 20)}...` : key)
        .style("font-family", "Roboto Condensed")
        .style("font-size", fontSize)
        .style("fill", "black")
    );
    pair = iter.next().value;
  }

  return centroids;
};
