/*****************************************
******************************************
      Setup SVG & Canvas elements
******************************************
******************************************/

const setupScatterplot = (
  width,
  height,
  margin
) => {

  var container = d3.select("#scatterplot")

  var svg = container.append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
    .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  var canvas = container.append("canvas")
    .attr("width", width * devicePixelRatio)
    .attr("height", height * devicePixelRatio)
    .style("width", width + "px")
    .style("height", height + "px")
    .style("margin-left", margin.left + 1 + "px") /* magic number: this seems to be visually correct when it is equal to the amount axis is translated by plus or minus... a little? */
    .style("margin-top", margin.top + "px")

  var ctx = canvas.node().getContext("2d");
      ctx.globalCompositeOperation = 'darken';
      ctx.globalAlpha = 0.15;
      ctx.lineWidth = 1.5;
      ctx.scale(devicePixelRatio, devicePixelRatio);

  return {
    svg,
    ctx,
  }

}

export default setupScatterplot;
