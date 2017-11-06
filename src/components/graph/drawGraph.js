import styles from "./graph.css";
import renderQueue from "../continuous/renderQueue";
import _ from "lodash";

const margin = {top: 20, right: 10, bottom: 30, left: 40},
    width = 960 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;

const x = d3.scaleLinear()
  .range([0, width]);

const y = d3.scaleLinear()
  .range([height, 0]);

/******************************************
*******************************************
          put svg & canvas in DOM
*******************************************
******************************************/

export const setupGraphElements = (data) => {

  var svg = d3.select("#graphAttachPoint").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .attr("class", `${styles.graphSVG}`)
    .append("g")
      .attr("transform", "translate(" + margin.left + " " + margin.top + ")");

  var chartArea = d3.select("#graphAttachPoint").append("div")
    .style("left", margin.left + "px")
    .style("top", margin.top + "px");

  var canvas = chartArea.append("canvas")
    .attr("width", width)
    .attr("height", height)
    .attr("class", `${styles.graphCanvas}`);

  var context = canvas.node().getContext("2d");

  return {
    svg,
    context
  }

}

/******************************************
*******************************************
        draw cells on the canvas
*******************************************
******************************************/

export const drawGraph = (data, context, expressionsCountsMap, color, ranges, metadata) => {

  /* clear canvas */
  context.clearRect(0, 0, width, height);

  let scale = null; /* it could be 'by expression' and that's a special case */

  if (
    color &&
    ranges[color].range /* set up a continuous scale */
  ) {
    scale = d3.scaleLinear()
      .domain([0, ranges[color].range.max])
      .range([1,0])
  }

  /* ! create a scale to map between expression values and colors, remove to somewhere else */
  // const expressionToColorScale = d3.scaleLinear()
  //   .domain([0, expressionsCountsMap.maxValue])
  //   .range([1,0])

  /* shuffle the data to overcome render order hiding cells */
  data = d3.shuffle(data);

  /* loop */
  data.forEach((p, i) => {
    context.beginPath();
    /* context.arc(x,y,r,sAngle,eAngle,counterclockwise); */
    context.arc(
      x(p[1]),            /* x */
      y(p[2]),            /* y */
      2,                  /* r */
      0,                  /* sAngle */
      2 * Math.PI         /* eAngle */
    );
    context.fill();

    if (i < 20) {
      // console.log(
      //   '0 to 1 scale', expressionToColorScale(expressionsCountsMap[p[0]]),
      //   'color', d3.interpolateViridis(expressionToColorScale(
      //     expressionsCountsMap[p[0]]
      //   ))
      // )
      // console.log('meta', _.find(metadata, {CellName: p[0]}), color, p)
    }
    // context.fillStyle = d3.interpolateViridis(expressionToColorScale(
    //   expressionsCountsMap[p[0]]
    // ));
    if (color) {
      context.fillStyle = d3.interpolateViridis(scale(
        _.find(metadata, {CellName: p[0]})[color] /* this.state.cells.metadata["23452345325"]["ERCC_reads"] = 20000 */
      ));
    }
  });
}
