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
  xscale,
  yAxis,
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
      const axes = drawAxes(
        this.state.svg,
        dimensions,
        xscale
      )
      drawParallelCoordinates(
        nextProps.metadata,
        dimensions,
        this.state.ctx,
        xscale,
        width,
        height,
        margin
      )
      this.setState({
        parallelExists: true,
        dimensions,
        axes,
      })
    }
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
