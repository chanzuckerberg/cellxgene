// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";

import HistogramBrush from "../brushableHistogram/";

@connect(state => {
  return {
    ranges: _.get(state.controls.world, "summary.obs", null),
    metadata: _.get(state.controls.world, "obsAnnotations", null),
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
    return () => {
      const { dispatch, ranges } = this.props;
      dispatch({
        type: "color by continuous metadata",
        colorAccessor: key,
        rangeMaxForColorAccessor: ranges[key].range.max
      });
    };
  }

  render() {
    const { ranges, obsAnnotations } = this.props;

    return (
      <div>
        {_.map(ranges, (value, key) => {
          const isColorField = key.includes("color") || key.includes("Color");
          if (value.range && key !== "name" && !isColorField) {
            return (
              <HistogramBrush
                key={key}
                field={key}
                fieldValues={obsAnnotations}
                ranges={value.range}
                handleColorAction={this.handleColorAction(key).bind(this)}
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
