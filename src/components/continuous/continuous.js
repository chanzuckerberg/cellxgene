import React from 'react';
import _ from "lodash";
import { connect } from "react-redux";

import styles from './parallelCoordinates.css';
import SectionHeader from "../framework/sectionHeader";

import renderQueue from "./renderQueue";

import setupParallelCoordinates from "./setupParallelCoordinates";
import drawAxes from "./drawAxes";
import drawLinesCanvas from "./drawLinesCanvas";

import {
  margin,
  width,
  height,
  innerHeight,
  color,
  createDimensions,
  types,
  processData,
  yAxis,
  d3_functor,
} from "./util";

@connect((state) => {

  const ranges = state.cells.cells && state.cells.cells.data.ranges ? state.cells.cells.data.ranges : null;
  const metadata = state.cells.cells && state.cells.cells.data.metadata ? state.cells.cells.data.metadata : null;

  return {
    ranges,
    metadata,
  }
})
class Continuous extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      svg: null,
      ctx: null,
      axes: null,
      parallelExists: false,
      dimensions: null,
    };
  }

  componentDidMount() {
    const {svg, ctx} = setupParallelCoordinates(
      width,
      height,
      margin
    );
    this.setState({svg, ctx})
  }

  componentWillReceiveProps(nextProps) {

    if (
      nextProps.ranges &&
      !this.state.parallelExists
    ) {
      const dimensions = createDimensions(nextProps.ranges);
      const xscale = d3.scalePoint()
        .domain(d3.range(dimensions.length))
        .range([0, width]);

      const {
        processedMetadata,
        processedDimensions
      } = processData(
        nextProps.metadata,
        dimensions
      )

      const _drawCellLines = this.drawCellLines(
        processedMetadata,
        processedDimensions,
        xscale
      );

      const axes = drawAxes(
        this.state.svg,
        this.state.ctx,
        processedDimensions,
        processedMetadata,
        xscale,
        height,
        width,
        _drawCellLines,
      )

      this.setState({
        parallelExists: true,
        dimensions,
        axes,
      })
    }
  }

  drawCellLines (metadata, dimensions, xscale) {
    const _drawCellLines = renderQueue(
      drawLinesCanvas(this.state.ctx, dimensions, xscale)
    ).rate(50);

    this.state.ctx.clearRect(0,0,width,height);
    this.state.ctx.globalAlpha = d3.min([0.85/Math.pow(metadata.length,0.3),1]);

    _drawCellLines(metadata);

    return _drawCellLines;

  }

  render() {

    return (
      <div id="parcoords_wrapper" style={{marginTop: 50}}>
        <SectionHeader text="Continuous Metadata"/>
        <div
          className={styles.parcoords}
          id="parcoords"
          style={{
            width:  width + margin.left + margin.right + "px",
            height: height + margin.top + margin.bottom + "px"
          }}></div>
      </div>
    )
  }
};

export default Continuous;
