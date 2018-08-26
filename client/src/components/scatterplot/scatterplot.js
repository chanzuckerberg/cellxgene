// jshint esversion: 6
// https://bl.ocks.org/Jverma/076377dd0125b1a508621441752735fc
// https://peterbeshai.com/scatterplot-in-d3-with-voronoi-interaction.html

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import _regl from "regl";
import * as d3 from "d3";

import _camera from "../../util/camera";

import setupScatterplot from "./setupScatterplot";
import styles from "./scatterplot.css";

import _drawPoints from "./drawPointsRegl";
import { scaleLinear } from "../../util/scaleLinear";

import { margin, width, height } from "./util";

@connect(state => {
  const { world } = state.controls2;
  const expressionX = world
    ? state.controls2.world.varDataByName(state.controls2.scatterplotXXaccessor)
    : null;
  const expressionY = world
    ? state.controls2.world.varDataByName(state.controls2.scatterplotYYaccessor)
    : null;

  return {
    world,

    colors: state.controls2.colors,
    colorAccessor: state.controls2.colorAccessor,
    colorScale: state.controls2.colorScale,

    // Accessors are var/gene names (strings)
    scatterplotXXaccessor: state.controls2.scatterplotXXaccessor,
    scatterplotYYaccessor: state.controls2.scatterplotYYaccessor,
    opacityForDeselectedCells: state.controls2.opacityForDeselectedCells,

    differential: state.differential,

    expressionX,
    expressionY,

    // updated whenever the crossfilter selection is updated
    selectionUpdate: _.get(state.controls2.world, "obsSelectionUpdateSeq", null)
  };
})
class Scatterplot extends React.Component {
  constructor(props) {
    super(props);
    this.count = 0;
    this.axes = false;
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
    let scales;
    const { expressionX, expressionY } = this.props;

    if (svg && expressionX && expressionY) {
      scales = Scatterplot.setupScales(expressionX, expressionY);
      this.drawAxesSVG(scales.xScale, scales.yScale, svg);
    }

    this.setState({
      svg,
      xScale: scales ? scales.xScale : null,
      yScale: scales ? scales.yScale : null
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

  componentDidUpdate(prevProps) {
    const {
      svg,
      xScale,
      yScale,
      regl,
      pointBuffer,
      colorBuffer,
      sizeBuffer
    } = this.state;
    const {
      world,
      scatterplotXXaccessor,
      scatterplotYYaccessor,
      expressionX,
      expressionY,
      colors
    } = this.props;

    if (
      world &&
      svg &&
      xScale &&
      yScale &&
      scatterplotXXaccessor &&
      scatterplotYYaccessor &&
      (scatterplotXXaccessor !== prevProps.scatterplotXXaccessor || // was CLU now FTH1 etc
      scatterplotYYaccessor !== prevProps.scatterplotYYaccessor || // was CLU now FTH1 etc
        !this.axes) // clicked off the tab and back again, rerender
    ) {
      this.drawAxesSVG(xScale, yScale, svg);
    }

    if (
      world &&
      regl &&
      pointBuffer &&
      colorBuffer &&
      sizeBuffer &&
      expressionX &&
      expressionY &&
      scatterplotXXaccessor &&
      scatterplotYYaccessor &&
      xScale &&
      yScale
    ) {
      const crossfilter = world.obsCrossfilter;
      const cellCount = expressionX.length;
      const positionsBuf = new Float32Array(2 * cellCount);
      const colorsBuf = new Float32Array(3 * cellCount);
      const sizesBuf = new Float32Array(cellCount);

      const glScaleX = scaleLinear([0, width], [-0.95, 0.95]);
      const glScaleY = scaleLinear([0, height], [-1, 1]);

      /*
        Construct Vectors
      */
      for (let i = 0; i < cellCount; i += 1) {
        positionsBuf[2 * i] = glScaleX(xScale(expressionX[i]));
        positionsBuf[2 * i + 1] = glScaleY(yScale(expressionY[i]));
      }

      for (let i = 0; i < cellCount; i += 1) {
        colorsBuf.set(colors.rgb[i], 3 * i);
      }

      crossfilter.fillByIsFiltered(sizesBuf, 4, 0.2);

      pointBuffer({ data: positionsBuf, dimension: 2 });
      colorBuffer({ data: colorsBuf, dimension: 3 });
      sizeBuffer({ data: sizesBuf, dimension: 1 });
      this.count = cellCount;
    }

    if (
      expressionX &&
      expressionY &&
      (this.props.scatterplotXXaccessor !== prevProps.scatterplotXXaccessor || // was CLU now FTH1 etc
        this.props.scatterplotYYaccessor !== prevProps.scatterplotYYaccessor)
    ) {
      const scales = Scatterplot.setupScales(expressionX, expressionY);
      this.setState(scales);
    }
  }

  static setupScales(expressionX, expressionY) {
    const xScale = d3
      .scaleLinear()
      .domain(d3.extent(expressionX))
      .range([0, width]);
    const yScale = d3
      .scaleLinear()
      .domain(d3.extent(expressionY))
      .range([height, 0]);

    return {
      xScale,
      yScale
    };
  }

  drawAxesSVG(xScale, yScale, svg) {
    svg.selectAll("*").remove();

    // the axes are much cleaner and easier now. No need to rotate and orient
    // the axis, just call axisBottom, axisLeft etc.
    const xAxis = d3.axisBottom().scale(xScale);

    const yAxis = d3.axisLeft().scale(yScale);

    // adding axes is also simpler now, just translate x-axis to (0,height)
    // and it's alread defined to be a bottom axis.
    svg
      .append("g")
      .attr("transform", "translate(0," + height + ")")
      .attr("class", "x axis")
      .call(xAxis);

    // y-axis is translated to (0,0)
    svg
      .append("g")
      .attr("transform", "translate(0,0)")
      .attr("class", "y axis")
      .call(yAxis);

    // adding label. For x-axis, it's at (10, 10), and for y-axis at (width, height-10).
    svg
      .append("text")
      .attr("x", 10)
      .attr("y", 10)
      .attr("class", "label")
      .text(this.props.scatterplotYYaccessor);

    svg
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
