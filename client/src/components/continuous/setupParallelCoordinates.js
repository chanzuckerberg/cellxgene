/*****************************************
******************************************
      Setup SVG & Canvas elements
******************************************
******************************************/
// jshint esversion: 6
const setupParallelCoordinates = (width, height, margin) => {
  var container = d3.select("#parcoords");

  var svg = container
    .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  var canvas = container
    .append("canvas")
    .attr("width", width * devicePixelRatio)
    .attr("height", height * devicePixelRatio)
    .style("width", width + "px")
    .style("height", height + "px")
    .style("margin-top", margin.top + "px")
    .style("margin-left", margin.left + "px");

  var ctx = canvas.node().getContext("2d");
  ctx.globalCompositeOperation = "darken";
  ctx.globalAlpha = 0.15;
  ctx.lineWidth = 1.5;
  ctx.scale(devicePixelRatio, devicePixelRatio);

  return {
    svg,
    ctx
  };
};

export default setupParallelCoordinates;
