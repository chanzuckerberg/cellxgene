// jshint esversion: 6
// https://bl.ocks.org/Jverma/076377dd0125b1a508621441752735fc
// https://peterbeshai.com/scatterplot-in-d3-with-voronoi-interaction.html

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";

import scatterplot from "./scatterplot";
import setupScatterplot from "./setupScatterplot";
import styles from "./scatterplot.css";

import mat4 from "gl-mat4";
import fit from "canvas-fit";
import _camera from "../../util/camera.js";
import _regl from "regl";
import _drawPoints from "./drawPointsRegl";
import { scaleLinear } from "../../util/scaleLinear";

import { margin, width, height, createDimensions } from "./util";

@connect(state => {
  const ranges =
    state.cells.cells && state.cells.cells.data.ranges
      ? state.cells.cells.data.ranges
      : null;
  const metadata =
    state.cells.cells && state.cells.cells.data.metadata
      ? state.cells.cells.data.metadata
      : null;
  const initializeRanges =
    state.initialize.data && state.initialize.data.data.ranges
      ? state.initialize.data.data.ranges
      : null;

  return {
    ranges,
    metadata,
    initializeRanges,
    colorAccessor: state.controls.colorAccessor,
    colorScale: state.controls.colorScale,
    scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
    scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
    opacityForDeselectedCells: state.controls.opacityForDeselectedCells,
    crossfilter: state.controls.crossfilter,
    differential: state.differential,
    expression: state.expression
  };
})
class Scatterplot extends React.Component {
  constructor(props) {
    super(props);
    this.count = 0;
    this.state = {
      svg: null,
      // ctx: null,
      axes: null,
      dimensions: null,
      xScale: null,
      yScale: null
    };
  }

  componentDidMount() {
    const { svg } = setupScatterplot(width, height, margin);
    this.setState({
      svg
    });

    const camera = _camera(this.reglCanvas, { scale: true, rotate: false });
    const regl = _regl(this.reglCanvas);

    const drawPoints = _drawPoints(regl);

    // preallocate buffers
    const pointBuffer = regl.buffer();
    const colorBuffer = regl.buffer();
    const sizeBuffer = regl.buffer();

    regl.frame(({ viewportWidth, viewportHeight }) => {
      regl.clear({
        depth: 1,
        color: [1, 1, 1, 1]
      });

      drawPoints({
        distance: camera.distance,
        color: colorBuffer,
        position: pointBuffer,
        size: sizeBuffer,
        count: this.count,
        view: camera.view(),
        scale: viewportHeight / viewportWidth
      });

      camera.tick();
    });

    this.setState({
      regl,
      sizeBuffer,
      pointBuffer,
      colorBuffer
    });
  }
  componentWillReceiveProps(nextProps) {
    this.maybeSetupScalesAndDrawAxes(nextProps);
  }
  componentDidUpdate(prevProps) {
    if (
      this.state.xScale &&
      this.state.yScale &&
      this.props.scatterplotXXaccessor &&
      this.props.scatterplotYYaccessor &&
      (this.props.scatterplotXXaccessor !== prevProps.scatterplotXXaccessor || // was CLU now FTH1 etc
        this.props.scatterplotYYaccessor !== prevProps.scatterplotYYaccessor)
    ) {
      this.drawAxesSVG(this.state.xScale, this.state.yScale);
    }

    if (
      this.state.regl &&
      this.state.pointBuffer &&
      this.state.colorBuffer &&
      this.state.sizeBuffer &&
      this.props.expression.data &&
      this.props.expression.data.genes &&
      this.props.scatterplotXXaccessor &&
      this.props.scatterplotYYaccessor &&
      this.state.xScale &&
      this.state.yScale
    ) {
      const crossfilter = this.props.crossfilter.cells;
      const data = this.props.expression.data;
      const cells = data.cells;
      const genes = data.genes;
      const cellCount = cells.length;
      const positions = new Float32Array(2 * cellCount);
      const colors = new Float32Array(3 * cellCount);
      const sizes = new Float32Array(cellCount);

      // d3.scaleLinear().domain([0, width]).range([-0.95, 0.95])
      const glScaleX = scaleLinear([0, width], [-0.95, 0.95]);

      // d3.scaleLinear().domain([0, height]).range([-1, 1])
      const glScaleY = scaleLinear([0, height], [-1, 1]);

      const geneXXaccessorIndex = genes.indexOf(
        this.props.scatterplotXXaccessor
      );
      const geneYYaccessorIndex = genes.indexOf(
        this.props.scatterplotYYaccessor
      );

      /*
        Construct Vectors
      */
      for (let i = 0; i < cellCount; i++) {
        const cell = cells[i];

        positions[2 * i] = glScaleX(
          this.state.xScale(cell.e[geneXXaccessorIndex])
        ); /* scale each point first to the window as we calculate extents separately below, so no need to repeat */
        positions[2 * i + 1] = glScaleY(
          this.state.yScale(cell.e[geneYYaccessorIndex])
        );
      }

      for (let i = 0; i < cellCount; i++) {
        const metadata = this.props.metadata[i];
        colors.set(metadata.__colorRGB__, 3 * i);
      }

      crossfilter.fillByIsFiltered(sizes, 4, 0.2);

      this.state.pointBuffer({ data: positions, dimension: 2 });
      this.state.colorBuffer({ data: colors, dimension: 3 });
      this.state.sizeBuffer({ data: sizes, dimension: 1 });
      this.count = cellCount;
    }
  }
  maybeSetupScalesAndDrawAxes(nextProps) {
    if (
      nextProps.expression &&
      nextProps.expression.data &&
      nextProps.scatterplotXXaccessor &&
      nextProps.scatterplotYYaccessor
    ) {
      const xScale = d3
        .scaleLinear()
        .domain(
          d3.extent(nextProps.expression.data.cells, (cell, i) => {
            return cell.e[
              nextProps.expression.data.genes.indexOf(
                nextProps.scatterplotXXaccessor
              )
            ];
          })
        )
        .range([0, width]);

      const yScale = d3
        .scaleLinear()
        .domain(
          d3.extent(nextProps.expression.data.cells, cell => {
            return cell.e[
              nextProps.expression.data.genes.indexOf(
                nextProps.scatterplotYYaccessor
              )
            ];
          })
        )
        .range([height, 0]);

      this.setState({
        xScale,
        yScale
      });
    }
  }
  drawAxesSVG(xScale, yScale) {
    this.state.svg.selectAll("*").remove();

    // the axes are much cleaner and easier now. No need to rotate and orient the axis, just call axisBottom, axisLeft etc.
    var xAxis = d3.axisBottom().scale(xScale);

    var yAxis = d3.axisLeft().scale(yScale);

    // adding axes is also simpler now, just translate x-axis to (0,height) and it's alread defined to be a bottom axis.
    this.state.svg
      .append("g")
      .attr("transform", "translate(0," + height + ")")
      .attr("class", "x axis")
      .call(xAxis);

    // y-axis is translated to (0,0)
    this.state.svg
      .append("g")
      .attr("transform", "translate(0,0)")
      .attr("class", "y axis")
      .call(yAxis);

    // adding label. For x-axis, it's at (10, 10), and for y-axis at (width, height-10).
    this.state.svg
      .append("text")
      .attr("x", 10)
      .attr("y", 10)
      .attr("class", "label")
      .text(this.props.scatterplotYYaccessor);

    this.state.svg
      .append("text")
      .attr("x", width)
      .attr("y", height - 10)
      .attr("text-anchor", "end")
      .attr("class", "label")
      .text(this.props.scatterplotXXaccessor);
  }

  render() {
    return (
      <div
        style={{
          backgroundColor: "white",
          paddingBottom: 20
        }}
        id="scatterplot_wrapper"
      >
        <div
          className={styles.scatterplot}
          id="scatterplot"
          style={{
            width: width + margin.left + margin.right + "px",
            height: height + margin.top + margin.bottom + "px"
          }}
        >
          <canvas
            width={width}
            height={height}
            style={{
              marginLeft: margin.left - 7,
              marginTop: margin.top
            }}
            ref={canvas => {
              this.reglCanvas = canvas;
            }}
          />
        </div>
      </div>
    );
  }
}

export default Scatterplot;

// <SectionHeader text="Continuous Metadata"/>
