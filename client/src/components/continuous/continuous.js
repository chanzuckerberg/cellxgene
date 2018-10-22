// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";

import HistogramBrush from "../brushableHistogram";

@connect(state => ({
  ranges: _.get(state.controls.world, "summary.obs", null),
  metadata: _.get(state.controls.world, "obsAnnotations", null),
  colorAccessor: state.controls.colorAccessor,
  colorScale: state.controls.colorScale,
  selectionUpdate: _.get(state.controls, "crossfilter.updateTime", null),
  schema: _.get(state.controls.world, "schema", null)
}))
class Continuous extends React.Component {
  constructor(props) {
    super(props);
    this.hasContinuous = false;
    this.continuousChecked = false;

    this.state = {};
  }

  componentDidUpdate() {}

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
    const { ranges, obsAnnotations, schema } = this.props;
    if (schema && !this.continuousChecked) {
      this.hasContinuous = _.some(
        schema.annotations.obs,
        d => d.type === "int32" || d.type === "float32"
      );
      this.continuousChecked = true; /* only do this once */
    }

    return (
      <div
        style={{
          padding: globals.leftSidebarSectionPadding
        }}
      >
        {this.hasContinuous ? (
          <p
            style={Object.assign({}, globals.leftSidebarSectionHeading, {
              marginTop: 40
            })}
          >
            Continuous metadata
          </p>
        ) : null}
        {_.map(ranges, (value, key) => {
          const isColorField = key.includes("color") || key.includes("Color");
          if (value.range && key !== "name" && !isColorField) {
            return (
              <HistogramBrush
                key={key}
                field={key}
                isObs
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
