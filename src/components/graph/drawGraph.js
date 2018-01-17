import styles from "./graph.css";
import renderQueue from "../continuous/renderQueue";
import _ from "lodash";
import * as globals from "../../globals";

/******************************************
*******************************************
          put svg & canvas in DOM
*******************************************
******************************************/

export const setupGraphElements = (
  handleBrushSelectAction,
  handleBrushDeselectAction
) => {

  var canvas = d3.select("#graphAttachPoint")
    .append("canvas")
    .attr("width", globals.graphWidth)
    .attr("height", globals.graphHeight)
    .attr("class", `${styles.graphCanvas}`);

  var ctx = canvas.node().getContext("2d");

  var svg = d3.select("#graphAttachPoint").append("svg")
    .attr("width", globals.graphWidth)
    .attr("height", globals.graphHeight)
    .attr("class", `${styles.graphSVG}`)
    // .append("g")
    //   .attr("transform", "translate(" + margin.left + " " + margin.top + ")");

  setupGraphBrush(
    svg,
    handleBrushSelectAction,
    handleBrushDeselectAction
  );

  return {
    svg,
    ctx
  }

}

/******************************************
*******************************************
        draw cells on the canvas
*******************************************
******************************************/

export const drawGraph = (
  context,
  expressionsCountsMap,
  colorAccessor,
  ranges,
  metadata,
  currentCellSelection,
  graphBrushSelection,
  colorScale,
  graphMap, /* tmp remove when structure exists on server */
  opacityForDeselectedCells,
) => {

  /* clear canvas */
  context.clearRect(0, 0, globals.graphWidth, globals.graphHeight);

  // let colorScale = null; /* it could be 'by expression' and that's a special case */

  /* ! create a scale to map between expression values and colors, remove to somewhere else */
  // const expressionToColorScale = d3.scaleLinear()
  //   .domain([0, expressionsCountsMap.maxValue])
  //   .range([1,0])

  const data = [];

  _.each(currentCellSelection, (cell, i) => {
    if (graphMap[cell["CellName"]]) { /* fails silently, sometimes this is undefined, in which case the graph array should be shorter than the cell array, check in reducer */
      data.push([
        cell["CellName"],
        graphMap[cell["CellName"]][0],
        graphMap[cell["CellName"]][1]
      ])
    }
  })

  /* shuffle the data to overcome render order hiding cells, & filter first */
  // data = d3.shuffle(data); /* make me a control */

  const _currentCellSelectionMap = _.keyBy(currentCellSelection, "CellName"); /* move me to the reducer */

  data.forEach((p, i) => {
    context.beginPath();
    /* context.arc(x,y,r,sAngle,eAngle,counterclockwise); */
    context.arc(
      globals.graphXScale(p[1]),            /* x */
      globals.graphYScale(p[2]),            /* y */
      2,                  /* r */
      0,                  /* sAngle */
      2 * Math.PI         /* eAngle */
    );

    context.fillStyle = _currentCellSelectionMap[p[0]]["__color__"]


    // if (colorAccessor && colorScale) {
    //   context.fillStyle = d3.interpolateViridis(colorScale(
    //     _currentCellSelectionMap[p[0]][colorAccessor]
    //   ));
    // } else {
    //   context.fillStyle = "rgb(0,0,0)";
    // }

    if (_currentCellSelectionMap[p[0]]["__selected__"]) {
      context.globalAlpha = 1;
    } else {
      context.globalAlpha = opacityForDeselectedCells;
    }

    // if (i < 20) {
      // console.log(
      //   '0 to 1 scale', expressionToColorScale(expressionsCountsMap[p[0]]),
      //   'color', d3.interpolateViridis(expressionToColorScale(
      //     expressionsCountsMap[p[0]]
      //   ))
      // )
      // console.log('meta', _.find(metadata, {CellName: p[0]}), color, p)
    // }
    // context.fillStyle = d3.interpolateViridis(expressionToColorScale(
    //   expressionsCountsMap[p[0]]
    // ));

    context.fill();
  });
}

const setupGraphBrush = (
  svg,
  handleBrushSelectAction,
  handleBrushDeselectAction
) => {
  svg.append("g")
   .call(
     d3.brush()
     .extent([
       [0, 0],
       [globals.graphWidth, globals.graphHeight]
     ])
     .on("brush", handleBrushSelectAction)
     .on("end", handleBrushDeselectAction)
   );
}
