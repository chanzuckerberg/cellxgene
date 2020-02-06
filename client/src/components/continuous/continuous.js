// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";
import { Button } from "@blueprintjs/core";

import HistogramBrush from "../brushableHistogram";

@connect(state => ({
  obsAnnotations: state.world?.obsAnnotations,
  colorAccessor: state.colors.colorAccessor,
  colorScale: state.colors.scale,
  schema: state.world?.schema
}))
class Continuous extends React.PureComponent {
  handleColorAction = key => {
    return () => {
      const { dispatch, obsAnnotations } = this.props;
      const summary = obsAnnotations.col(key).summarize();
      dispatch({
        type: "color by continuous metadata",
        colorAccessor: key,
        rangeForColorAccessor: summary
      });
    };
  };

  renderIsStillLoading(zebra, key) {
    return (
      <div
        key={key}
        style={{
          padding: globals.leftSidebarSectionPadding,
          backgroundColor: zebra % 2 === 0 ? globals.lightestGrey : "white"
        }}
      >
        <div
          style={{
            display: "flex",
            justifyContent: "space-between",
            justifyItems: "center",
            alignItems: "center"
          }}
        >
          <div style={{ minWidth: 30 }}></div>
          <div style={{ display: "flex", alignSelf: "center" }}>
            <span style={{ fontStyle: "italic" }}>{key}</span>
          </div>
          <div
            style={{
              display: "flex",
              justifyContent: "flex-end"
            }}
          >
            <Button minimal loading intent="primary" />
          </div>
        </div>
      </div>
    );
  }

  render() {
    const { obsAnnotations, schema } = this.props;

    const obsIndex = schema.annotations.obs.index;
    const allContinuousNames = schema.annotations.obs.columns
      .filter(col => col.type === "int32" || col.type === "float32")
      .filter(col => col.name != obsIndex)
      .map(col => col.name);

    /* initial value for iterator to simulate index, ranges is an object */
    let zebra = 0;

    return (
      <div>
        {allContinuousNames.map(key => {
          if (!obsAnnotations.hasCol(key)) {
            // still loading!
            zebra += 1;
            return this.renderIsStillLoading(zebra, key);
          } else {
            // data loaded and available
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
                  handleColorAction={this.handleColorAction(key)}
                />
              );
            }
          }
        })}
      </div>
    );
  }
}

export default Continuous;
