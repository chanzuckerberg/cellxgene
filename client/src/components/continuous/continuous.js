// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";

import styles from "./parallelCoordinates.css";
import SectionHeader from "../framework/sectionHeader";

import setupParallelCoordinates from "./setupParallelCoordinates";
import drawAxes from "./drawAxes";
import drawLinesCanvas from "./drawLinesCanvas";

import HistogramBrush from "./histogramBrush";

import { margin, width, height, createDimensions } from "./util";

@connect(state => {
  const metadata = _.get(state.controls2.world, "obsAnnotations", null);
  const ranges = _.get(state.controls2.world, "summary.obs", null);

  return {
    ranges,
    metadata,
    colorAccessor: state.controls2.colorAccessor,
    colorScale: state.controls2.colorScale,
    selectionUpdate: _.get(state.controls2, "obsCrossfilter.updateTime", null)
  };
})
class Continuous extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      svg: null,
      ctx: null,
      axes: null,
      dimensions: null
    };
  }

  handleBrushAction(selection) {
    this.props.dispatch({
      type: "continuous selection using parallel coords brushing",
      data: selection
    });
  }

  handleColorAction(key) {
    this.props.dispatch({
      type: "color by continuous metadata",
      colorAccessor: key,
      rangeMaxForColorAccessor: this.props.ranges[key].range.max
    });
  }

  render() {
    return (
      <div>
        {_.map(this.props.ranges, (value, key) => {
          const isColorField = key.includes("color") || key.includes("Color");
          if (value.range && key !== "name" && !isColorField) {
            return (
              <HistogramBrush
                key={key}
                metadataField={key}
                ranges={value.range}
              />
            );
          }
        })}
      </div>
    );
  }
}

export default Continuous;

// <SectionHeader text="Continuous Metadata"/>
