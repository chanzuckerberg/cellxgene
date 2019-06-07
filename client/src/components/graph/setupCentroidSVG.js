import * as d3 from "d3";
import styles from "./graph.css";

export default (responsive, graphPaddingRight, xy) => {
  const svg = d3
    .select("#graphAttachPoint")
    .append("svg")
    .attr("id", "centroid")
    .attr("data-testid", "centroid-overlay")
    .attr("width", responsive.width - graphPaddingRight)
    .attr("height", responsive.height)
    .attr("class", `${styles.graphSVG}`)
    .style("z-index", 998);
  //  TODO: Create own styles, ask Colin for an explanation on the css
  // For now I'm going to put centroid z-index at 998 and lasso on 999

  const offset = 5;

  svg
    .append("circle")
    .attr("cx", xy[0] + offset)
    .attr("cy", xy[1] + offset)
    .attr("r", 10)
    .style("fill", "rgba(244, 66, 66, .55)")
    .style("stroke", "rgba(0, 0, 0, .4");

  return svg;
};
