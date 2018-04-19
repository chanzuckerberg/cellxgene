// https://bl.ocks.org/Jverma/076377dd0125b1a508621441752735fc
// https://peterbeshai.com/scatterplot-in-d3-with-voronoi-interaction.html

import React from 'react';
import _ from "lodash";
import { connect } from "react-redux";

import scatterplot from "./scatterplot";
import setupScatterplot from "./setupScatterplot";
import styles from './scatterplot.css';
import drawScatterplotCanvas from "./drawScatterplotCanvas";

import mat4 from 'gl-mat4';
import fit from 'canvas-fit';
import _camera from '../../util/camera.js'
import _regl from 'regl'
import _drawPoints from './drawPointsRegl'

import {
  scaleRGB
} from "../../util/scaleRGB";

import {
  margin,
  width,
  height,
  createDimensions,
} from "./util";

@connect((state) => {

  const ranges = state.cells.cells && state.cells.cells.data.ranges ? state.cells.cells.data.ranges : null;
  const metadata = state.cells.cells && state.cells.cells.data.metadata ? state.cells.cells.data.metadata : null;
  const initializeRanges = state.initialize.data && state.initialize.data.data.ranges ? state.initialize.data.data.ranges : null;

  return {
    ranges,
    metadata,
    initializeRanges,
    colorAccessor: state.controls.colorAccessor,
    colorScale: state.controls.colorScale,
    currentCellSelection: state.controls.currentCellSelection,
    scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
    scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
    opacityForDeselectedCells: state.controls.opacityForDeselectedCells,
    differential: state.differential,
    expression: state.expression,
  }
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
      yScale: null,
    };
  }

  componentDidMount() {
    const {
      svg
      // ctx
    } = setupScatterplot(
      width,
      height,
      margin
    );
    this.setState({
      svg
      // ctx
    })

    const camera = _camera(this.reglCanvas, {scale: true, rotate: false});
    const regl = _regl(this.reglCanvas)

    const drawPoints = _drawPoints(regl)

    // preallocate buffers
    const pointBuffer = regl.buffer()
    const colorBuffer = regl.buffer()
    const sizeBuffer = regl.buffer()


    regl.frame(({viewportWidth, viewportHeight}) => {

      regl.clear({
        depth: 1,
        color: [1, 1, 1, 1]
      })

      drawPoints({
        distance: camera.distance,
        color: colorBuffer,
        position: pointBuffer,
        size: sizeBuffer,
        count: this.count,
        view: camera.view(),
        scale: viewportHeight / viewportWidth
      })

      camera.tick()
    })

    this.setState({
      regl,
      sizeBuffer,
      pointBuffer,
      colorBuffer,
    })

  }
  componentWillReceiveProps(nextProps) {

    this.maybeSetupScalesAndDrawAxes(nextProps);

  }
  componentDidUpdate(prevProps) {
    if (
      (this.state.xScale && this.state.yScale) &&
      (this.props.scatterplotXXaccessor && this.props.scatterplotYYaccessor) &&
      this.props.scatterplotXXaccessor !== prevProps.scatterplotXXaccessor || // was CLU now FTH1 etc
      this.props.scatterplotYYaccessor !== prevProps.scatterplotYYaccessor
    ) {
      this.drawAxesSVG(this.state.xScale, this.state.yScale);
    }

    // if (
    //   this.state.xScale &&
    //   this.state.yScale
    // ) {
    //   drawScatterplotCanvas(
    //     this.state.ctx,
    //     this.state.xScale,
    //     this.state.yScale,
    //     this.props.currentCellSelection,
    //     this.props.opacityForDeselectedCells,
    //     this.props.expression,
    //     this.props.scatterplotXXaccessor,
    //     this.props.scatterplotYYaccessor,
    //   )
    // }

    if (
      this.state.regl &&
      this.state.pointBuffer &&
      this.state.colorBuffer &&
      this.state.sizeBuffer &&
      this.props.currentCellSelection &&
      this.props.expression.data &&
      this.props.expression.data.genes &&
      this.props.scatterplotXXaccessor &&
      this.props.scatterplotYYaccessor &&
      this.state.xScale &&
      this.state.yScale
    ) {
      const _currentCellSelectionMap = _.keyBy(this.props.currentCellSelection, "CellName"); /* move me to the reducer */

      const positions = [];
      const colors = [];
      const sizes = [];

      const glScaleX = d3.scaleLinear()
        .domain([0, width])
        .range([-.95, .95]) /* padding */

      const glScaleY = d3.scaleLinear()
        .domain([0, height])
        .range([-1, 1])


      /*
        Construct Vectors
      */
      _.each(this.props.expression.data.cells, (cell, i) => {
        /*
          this if is necessary until we are no longer getting expression for all cells, but only for 'world'
          ...which will mean refetching when we regraph, or 'go back up to all cells'
        */
        if (_currentCellSelectionMap[cell.cellname]) { /* fails silently, sometimes this is undefined, in which case the graph array should be shorter than the cell array, check in reducer */
          positions.push([
            glScaleX(this.state.xScale(cell.e[this.props.expression.data.genes.indexOf(this.props.scatterplotXXaccessor)])), /* scale each point first to the window as we calculate extents separately below, so no need to repeat */
            glScaleY(this.state.yScale(cell.e[this.props.expression.data.genes.indexOf(this.props.scatterplotYYaccessor)]))
          ])

          let c = _currentCellSelectionMap[cell.cellname]["__color__"];

          if (c[0] !== "#") {
            const _c = c.replace(/[^\d,.]/g, '').split(',');
            colors.push([
              scaleRGB(+_c[0]),
              scaleRGB(+_c[1]),
              scaleRGB(+_c[2])
            ])
          } else {
            var parsedHex = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(c);
            colors.push([
              scaleRGB(parseInt(parsedHex[1], 16)),
              scaleRGB(parseInt(parsedHex[2], 16)),
              scaleRGB(parseInt(parsedHex[3], 16))
            ]);
          }
          sizes.push(_currentCellSelectionMap[cell.cellname]["__selected__"] ? 4 : .2) /* make this a function of the number of total cells, including regraph */
        }
      })

      this.state.pointBuffer(positions)
      this.state.colorBuffer(colors)
      this.state.sizeBuffer(sizes)
      // this.state.pointBuffer(generatePoints(count))
      // this.state.colorBuffer(generateColors(count))
      // this.state.sizeBuffer(generateSizes(count))
      this.count = positions.length
    }
  }
  maybeSetupScalesAndDrawAxes(nextProps) {
    if (
      nextProps.expression &&
      nextProps.expression.data &&
      nextProps.scatterplotXXaccessor &&
      nextProps.scatterplotYYaccessor
    ) {
      const xScale = d3.scaleLinear()
        .domain(d3.extent(nextProps.expression.data.cells, (cell, i) => {
          return cell.e[nextProps.expression.data.genes.indexOf(nextProps.scatterplotXXaccessor)]
        }))
        .range([0, width])

      const yScale = d3.scaleLinear()
        .domain(d3.extent(nextProps.expression.data.cells, (cell) => {
          return cell.e[nextProps.expression.data.genes.indexOf(nextProps.scatterplotYYaccessor)]
        }))
        .range([height, 0])

      this.setState({
        xScale,
        yScale
      })
    }
  }
  drawAxesSVG(xScale, yScale) {

    this.state.svg.selectAll("*").remove();

    // the axes are much cleaner and easier now. No need to rotate and orient the axis, just call axisBottom, axisLeft etc.
    var xAxis = d3.axisBottom()
      .scale(xScale);

    var yAxis = d3.axisLeft()
      .scale(yScale);

    // adding axes is also simpler now, just translate x-axis to (0,height) and it's alread defined to be a bottom axis.
    this.state.svg.append('g')
      .attr('transform', 'translate(0,' + height + ')')
      .attr('class', 'x axis')
      .call(xAxis);

    // y-axis is translated to (0,0)
    this.state.svg.append('g')
      .attr('transform', 'translate(0,0)')
      .attr('class', 'y axis')
      .call(yAxis);

    // adding label. For x-axis, it's at (10, 10), and for y-axis at (width, height-10).
    this.state.svg.append('text')
      .attr('x', 10)
      .attr('y', 10)
      .attr('class', 'label')
      .text(this.props.scatterplotYYaccessor);

    this.state.svg.append('text')
      .attr('x', width)
      .attr('y', height - 10)
      .attr('text-anchor', 'end')
      .attr('class', 'label')
      .text(this.props.scatterplotXXaccessor);

  }

  render() {
    return (
      <div
        style={{
          backgroundColor: "white",
          borderRadius: 3,
          marginTop: 15,
          paddingBottom: 20,
          boxShadow: "3px 4px 13px 0px rgba(201,201,201,1)",
        }}
        id="scatterplot_wrapper">
        <div
          className={styles.scatterplot}
          id="scatterplot"
          style={{
            width:  width + margin.left + margin.right + "px",
            height: height + margin.top + margin.bottom + "px",
          }}
        >
          <canvas
            width={width}
            height={height}
            style={{
              marginLeft: margin.left - 5,
              marginTop: margin.top
            }}
            ref={(canvas) => { this.reglCanvas = canvas}}/>
        </div>
      </div>
    )
  }
};

export default Scatterplot;


// <SectionHeader text="Continuous Metadata"/>
