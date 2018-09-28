// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";

import HistogramBrush from "./histogramBrush";

@connect(state => {
  const metadata = _.get(state.controls.world, "obsAnnotations", null);
  const ranges = _.get(state.controls.world, "summary.obs", null);

  return {
    ranges,
    metadata,
    colorAccessor: state.controls.colorAccessor,
    colorScale: state.controls.colorScale,
    selectionUpdate: _.get(state.controls, "crossfilter.updateTime", null)
  };
})
class Continuous extends React.Component {
  handleBrushAction(selection) {
    const { dispatch } = this.props;
    dispatch({
      type: "continuous selection using parallel coords brushing",
      data: selection
    });
  }

  handleColorAction(key) {
    const { dispatch, ranges } = this.props;
    dispatch({
      type: "color by continuous metadata",
      colorAccessor: key,
      rangeMaxForColorAccessor: ranges[key].range.max
    });
  }

  render() {
    const { ranges } = this.props;

    return (
      <div>
        {_.map(ranges, (value, key) => {
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
          return null;
        })}
      </div>
    );
  }
}

export default Continuous;
