// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";

import HistogramBrush from "../brushableHistogram";

@connect(state => ({
  obsAnnotations: _.get(state.controls.world, "obsAnnotations", null),
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
      const { dispatch, obsAnnotations } = this.props;
      const summary = obsAnnotations.col(key).summarize();
      dispatch({
        type: "color by continuous metadata",
        colorAccessor: key,
        rangeMaxForColorAccessor: summary.max
      });
    };
  }

  render() {
    const { obsAnnotations, schema } = this.props;
    if (schema && !this.continuousChecked) {
      this.hasContinuous = _.some(
        schema.annotations.obs,
        d => d.type === "int32" || d.type === "float32"
      );
      this.continuousChecked = true; /* only do this once */
    }

    /* initial value for iterator to simulate index, ranges is an object */
    let zebra = -1;

    return (
      <div>
        {this.hasContinuous ? (
          <p
            style={Object.assign({}, globals.leftSidebarSectionHeading, {
              marginTop: 40,
              paddingLeft: globals.leftSidebarSectionPadding
            })}
          >
            Continuous metadata
          </p>
        ) : null}
        {obsAnnotations
          ? _.map(obsAnnotations.colIndex.keys(), key => {
              const summary = obsAnnotations.col(key).summarize();
              const isColorField =
                key.includes("color") || key.includes("Color");
              const nonFiniteExtent =
                summary.min === undefined || summary.max === undefined;
              zebra += 1;
              if (
                !summary.categorical &&
                key !== "name" &&
                !isColorField &&
                !nonFiniteExtent
              ) {
                return (
                  <HistogramBrush
                    key={key}
                    field={key}
                    isObs
                    zebra={zebra % 2 === 0}
                    ranges={summary}
                    handleColorAction={this.handleColorAction(key).bind(this)}
                  />
                );
              }
              return null;
            })
          : null}
      </div>
    );
  }
}

export default Continuous;
