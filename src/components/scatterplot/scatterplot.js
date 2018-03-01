// https://bl.ocks.org/Jverma/076377dd0125b1a508621441752735fc
// https://peterbeshai.com/scatterplot-in-d3-with-voronoi-interaction.html

import React from 'react';
import _ from "lodash";
import { connect } from "react-redux";

import scatterplot from "./scatterplot";
import setupScatterplot from "./setupScatterplot";
import styles from './scatterplot.css';
import drawScatterplotCanvas from "./drawScatterplotCanvas";

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
    this.state = {
      svg: null,
      ctx: null,
      axes: null,
      dimensions: null,
      xScale: null,
      yScale: null,
    };
  }

  componentDidMount() {
    const {svg, ctx} = setupScatterplot(
      width,
      height,
      margin
    );
    this.setState({svg, ctx})
  }
  componentWillReceiveProps(nextProps) {
    this.maybeSetupScalesAndDrawAxes(nextProps);
  }
  componentDidUpdate(prevProps) {
    if (
      (this.props.scatterplotXXaccessor && this.props.scatterplotYYaccessor) &&
      this.props.scatterplotXXaccessor !== prevProps.scatterplotXXaccessor || // was CLU now FTH1 etc
      this.props.scatterplotYYaccessor !== prevProps.scatterplotYYaccessor
    ) {
      this.drawAxesSVG(this.state.xScale, this.state.yScale);
    }

    if (
      this.state.xScale &&
      this.state.yScale
    ) {
      /* clear canvas */
      this.state.ctx.clearRect(0, 0, width, height);
      drawScatterplotCanvas(
        this.state.ctx,
        this.state.xScale,
        this.state.yScale,
        this.props.currentCellSelection,
        this.props.opacityForDeselectedCells,
        this.props.expression,
        this.props.scatterplotXXaccessor,
        this.props.scatterplotYYaccessor,
      )
    }
  }
  maybeSetupScalesAndDrawAxes(nextProps) {
    if (
      nextProps.expression &&
      nextProps.scatterplotXXaccessor &&
      nextProps.scatterplotYYaccessor
    ) {
      console.log('expression in scatterplot: ', nextProps)
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
        ></div>
      </div>
    )
  }
};

export default Scatterplot;


// <SectionHeader text="Continuous Metadata"/>
