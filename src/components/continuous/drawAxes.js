import styles from './parallelCoordinates.css';
import {
  yAxis,
  brushstart,
} from "./util";


const drawAxes = (
  svg,
  ctx,
  dimensions,
  metadata,
  xscale,
  height,
  width,
  handleBrushAction,
) => {

  /*****************************************
  ******************************************
  Handles a brush event, toggling the display of foreground lines.
  ******************************************
  ******************************************/

  function brush () {
    // _drawCellLines.invalidate(); /* this should be moved up */

    var actives = [];
    svg.selectAll(".parcoords_axis .parcoords_brush")
      .filter(function (d) {
        return d3.brushSelection(this);
      })
      .each(function(d) {
        actives.push({
          dimension: d,
          extent: d3.brushSelection(this)
        });
      });

    var selected = metadata.filter(function(d) {
      /* this is iterating over the enter dataset */
      if (actives.every(function(active) {
          var dim = active.dimension;
          // test if point is within extents for each active brush
          return dim.type.within(d[dim.key], active.extent, dim);
        })) {
        return true;
      }
    });

    /* send the selected cells to redux */
    handleBrushAction(selected)

    /* Reset canvas */
    // ctx.clearRect(0, 0, width, height);
    // ctx.globalAlpha = d3.min([0.85 / Math.pow(selected.length, 0.3), 1]);

    /* pass the result of the filter above to a fresh canvas, using RAF */
    // _drawCellLines(selected);

    /* text beneath */
    // output.text(d3.tsvFormat(selected.slice(0,22)));
  }

  var axes = svg.selectAll(".parcoords_axis")
      .data(dimensions)
    .enter().append("g")
      .attr("class", `${styles.axis} parcoords_axis`)
      .attr("transform", (d,i) => { return "translate(" + xscale(i) + ")"; });

  axes.append("g")
      .each(function(d) {
        var renderAxis = "axis" in d
          ? d.axis.scale(d.scale)  // custom axis
          : yAxis.scale(d.scale);  // default axis
        d3.select(this).call(renderAxis);
      })
    .append("text")
      .attr("class", styles.title)
      .attr("text-anchor", "start")
      .text(function(d) { return "description" in d ? d.description : d.key; });

  // Add and store a brush for each axis.
  axes.append("g")
      .attr("class", `${styles.brush} parcoords_brush`)
      .each(function(d) {
        d3.select(this).call(
          d.brush = d3.brushY()
            .extent([[-10,0], [10, height]])
            .on("start", brushstart)
            .on("brush", brush)
            .on("end", brush)
        )
      })
    .selectAll("rect")
      .attr("x", -8)
      .attr("width", 16);

  return axes;

}

export default drawAxes;
