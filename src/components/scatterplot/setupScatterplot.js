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


  return {
    svg,
  }

}

export default setupScatterplot;
