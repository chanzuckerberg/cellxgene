import React from 'react';
import _ from "lodash";
import { connect } from "react-redux";

import styles from './parallelCoordinates.css';
import SectionHeader from "../framework/sectionHeader";

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

  const initializeRanges = state.initialize.data && state.initialize.data.data.ranges ? state.initialize.data.data.ranges : null;
  const initializeMetadata = state.initialize.data && state.initialize.data.data.metadata ? state.initialize.data.data.metadata : null;

  return {
    ranges,
    metadata,
    initializeRanges,
    initializeMetadata,
    color: state.controls.color,
    continuousSelection: state.controls.continuousSelection
  }
})
class Continuous extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      svg: null,
      ctx: null,
      axes: null,
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
    this.maybeDrawLines(nextProps);
    this.maybeDrawAxes(nextProps);
  }

  maybeDrawAxes(nextProps) {
    if (
      !this.state.axes &&
      nextProps.initializeRanges /* we assume if this is present, everything else is too */
    ) {
      const dimensions = createDimensions(nextProps.initializeRanges);
      const xscale = d3.scalePoint()
        .domain(d3.range(dimensions.length))
        .range([0, width]);

      const {
        processedMetadata,
        processedDimensions
      } = processData(
        nextProps.initializeMetadata,
        dimensions
      )

      const axes = drawAxes(
        this.state.svg,
        this.state.ctx,
        processedDimensions,
        processedMetadata,
        xscale,
        height,
        width,
        this.handleBrushAction.bind(this),
      );

      this.setState({
        axes,
        xscale,
        processedMetadata,
        processedDimensions,
      })

    }
  }

  maybeDrawLines(nextProps) {
    if (
      nextProps.ranges &&
      nextProps.continuousSelection &&
      this.state.axes /* this may cause things to never render :/ forceUpdate? */
    ) {
      if (this.state._drawLinesCanvas) {
        this.state._drawLinesCanvas.invalidate();
        this.state.ctx.clearRect(0,0,width,height);
        this.state.ctx.globalAlpha = d3.min([
          0.85 / Math.pow(this.props.initializeMetadata.length, 0.3),
          1
        ]);
      }

      console.log(        _.filter(this.state.processedMetadata, (d) => {
                return nextProps.continuousSelection.indexOf(d.CellName) > -1
              }))

      const _drawLinesCanvas = drawLinesCanvas(
        _.filter(this.state.processedMetadata, (d) => {
          return nextProps.continuousSelection.indexOf(d.CellName) > -1
        }),
        this.state.processedDimensions,
        this.state.xscale,
        this.state.ctx,
        width,
        height
      );

      this.setState({
        // dimensions,
        _drawLinesCanvas,
      })
    }
  }

  handleBrushAction (selection) {
    this.props.dispatch({
      type: "continuous selection using parallel coords brushing",
      data: _.map(selection, (d) => { return d.CellName })
    })
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
