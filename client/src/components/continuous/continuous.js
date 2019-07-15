// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";

import HistogramBrush from "../brushableHistogram";

@connect(state => ({
  obsAnnotations: state.world?.obsAnnotations,
  colorAccessor: state.colors.colorAccessor,
  colorScale: state.colors.scale,
  schema: state.world?.schema
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
        rangeForColorAccessor: summary
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
    let zebra = 0;

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
              const isColorField =
                key.includes("color") || key.includes("Color");
              if (key === schema.annotations.obs.index || isColorField)
                return null;

              const summary = obsAnnotations.col(key).summarize();
              const nonFiniteExtent =
                summary.min === undefined ||
                summary.max === undefined ||
                Number.isNaN(summary.min) ||
                Number.isNaN(summary.max);
              if (!summary.categorical && !nonFiniteExtent) {
                zebra += 1;
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
