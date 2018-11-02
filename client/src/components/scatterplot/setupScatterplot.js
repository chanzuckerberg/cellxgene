// jshint esversion: 6
/*****************************************
******************************************
      Setup SVG & Canvas elements
******************************************
******************************************/

import * as d3 from "d3";

const setupScatterplot = (width, height, margin) => {
  const container = d3.select("#scatterplot");

  const svg = container
    .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);

  return {
    svg
  };
};

export default setupScatterplot;
