import * as d3 from "d3";
import styles from "./graph.css";

export default (
  responsive,
  graphPaddingRight
) => {
  const svg = d3
  .select("#graphAttachPoint")
  .append("svg")
  .attr("id", "centroid")
  .attr("data-testid", "centroid-overlay")
  .attr("width", responsive.width - graphPaddingRight)
  .attr("height", responsive.height)
  .attr("class", `${styles.graphSVG}`);
//  TODO: Create own styles, ask Colin for an explanation on the css
  return svg;
}