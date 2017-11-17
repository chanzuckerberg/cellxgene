import styles from "./graph.css";
import renderQueue from "../continuous/renderQueue";
import _ from "lodash";

const margin = {top: 20, right: 10, bottom: 30, left: 40},
    width = 960 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;

const x = d3.scaleLinear()
  .domain([0, 1]) /* while this is the default for d3, our data is normalized so better to be explicit */
  .range([0, width]);

const y = d3.scaleLinear()
  .domain([0, 1]) /* while this is the default for d3, our data is normalized so better to be explicit */
  .range([height, 0]);

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
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .attr("class", `${styles.graphCanvas}`);

  var ctx = canvas.node().getContext("2d");

  var svg = d3.select("#graphAttachPoint").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
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
  data,
  context,
  expressionsCountsMap,
  colorAccessor,
  ranges,
  metadata,
  continuousSelection,
  graphBrushSelection,
  colorScale,
) => {

  /* clear canvas */
  context.clearRect(0, 0, width, height);

  // let colorScale = null; /* it could be 'by expression' and that's a special case */

  /* ! create a scale to map between expression values and colors, remove to somewhere else */
  // const expressionToColorScale = d3.scaleLinear()
  //   .domain([0, expressionsCountsMap.maxValue])
  //   .range([1,0])

  if (continuousSelection) {
    data = _.filter(data, (d) => {return continuousSelection.indexOf(d[0]) > -1 })
  }

  /* shuffle the data to overcome render order hiding cells, & filter first */
  data = d3.shuffle(data); /* make me a control */

  const _metadata = _.keyBy(metadata, "CellName"); /* move me to the reducer */

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

    let pointIsInsideBrushBounds;

    if (graphBrushSelection) {

      pointIsInsideBrushBounds = (
        x(p[1]) >= graphBrushSelection.northwestX &&
        x(p[1]) <= graphBrushSelection.southeastX &&
        y(p[2]) >= graphBrushSelection.northwestY &&
        y(p[2]) <= graphBrushSelection.southeastY
      );

      // if (i === 5) {
      //   console.log(
      //     pointIsInsideBrushBounds,
      //     x(p[1]),
      //     y(p[2]),
      //     graphBrushSelection,
      //   );
      // }
    }

    context.globalAlpha = !graphBrushSelection || pointIsInsideBrushBounds ? 1 : .2;

    if (colorAccessor && colorScale) {
      context.fillStyle = d3.interpolateViridis(colorScale(
        _metadata[p[0]][colorAccessor]
        // _.find(metadata, {CellName: p[0]})[color] /* this.state.cells.metadata["23452345325"]["ERCC_reads"] = 20000 would be much faster as key value lookup */
      ));
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
       [width, height]
     ])
     .on("brush", handleBrushSelectAction)
     .on("end", handleBrushDeselectAction)
   );
}
