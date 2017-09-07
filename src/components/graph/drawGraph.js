var margin = {top: 20, right: 10, bottom: 30, left: 40},
    width = 960 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;

var svg = d3.select("body").append("svg")
  .attr("width", width + margin.left + margin.right)
  .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform", "translate(" + margin.left + " " + margin.top + ")");

var factory = d3.geom.quadtree()
  .extent([
    [0, 0],
    [width, height]
  ]);

var x = d3.scale.linear()
  .range([0, width]);

var y = d3.scale.linear()
  .range([height, 0]);

var xAxis = d3.svg.axis()
  .scale(x)
  .orient("bottom");

var yAxis = d3.svg.axis()
  .scale(y)
  .orient("left");

var xg = svg.append("g")
  .attr("class", "x axis")
  .attr("transform", "translate(0," + height + ")");

var yg = svg.append("g")
  .attr("class", "y axis");

var chartArea = d3.select("body").append("div")
  .style("left", margin.left + "px")
  .style("top", margin.top + "px");

var canvas = chartArea.append("canvas")
  .attr("width", width)
  .attr("height", height);

var context = canvas.node().getContext("2d");

context.fillStyle = "#f0f";

// Layer on top of canvas, example of selection details
var highlight = chartArea.append("svg")
  .attr("width", width)
  .attr("height", height)
    .append("circle")
      .attr("r", 7)
      .classed("hidden", true);

redraw();

function redraw() {

  // Randomize the scale
  var scale = 1 + Math.floor(Math.random() * 10);

  // Redraw axes
  x.domain([0, scale]);
  y.domain([0, scale]);
  xg.call(xAxis);
  yg.call(yAxis);

  var points = randomPoints(scale);

  var tree = factory(points);

  // Update canvas
  context.clearRect(0, 0, width, height);

  points.forEach(function(p,i){

    context.beginPath();
    context.arc(x(p[0]), y(p[1]), 5, 0, 2 * Math.PI);
    context.fill();

  });

  canvas.on("mousemove",function(){

    var mouse = d3.mouse(this),
        closest = tree.find([x.invert(mouse[0]), y.invert(mouse[1])]);

    highlight.attr("cx", x(closest[0]))
      .attr("cy", y(closest[1]));

  });

  canvas.on("mouseover",function(){
    highlight.classed("hidden", false);
  });

  canvas.on("mouseout",function(){
    highlight.classed("hidden", true);
  });

}

function randomPoints(scale) {

  // Get points
  return d3.range(1000).map(function(d){

    return [
      Math.random() * scale,
      Math.random() * scale
    ];

  });

}
