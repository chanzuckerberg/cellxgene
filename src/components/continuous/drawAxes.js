import styles from './parallelCoordinates.css';

const drawAxes = (
  svg,
  dimensions,
  xscale
) => {

  var axes = svg.selectAll(".parcoords_axis")
      .data(dimensions)
    .enter().append("g")
      .attr("class", `${styles.axis} parcoords_axis`)
      .attr("transform", (d,i) => { return "translate(" + xscale(i) + ")"; });

  return axes;

}

export default drawAxes;
