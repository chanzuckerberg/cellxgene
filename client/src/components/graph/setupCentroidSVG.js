import * as d3 from "d3";

export default (
  responsive,
  graphPaddingRight
) => {
  const svg = d3
  .select("#graphAttachPoint")
  .append("svg")
  .attr("data-testid", "centroid-overlay")
  .attr("width", responsive.width - graphPaddingRight)
  .attr("height", responsive.height);

  return svg;
}