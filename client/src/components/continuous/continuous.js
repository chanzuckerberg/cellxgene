// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import { connect } from "react-redux";
import { Button } from "@blueprintjs/core";

import * as globals from "../../globals";
import HistogramBrush from "../brushableHistogram";

@connect((state) => ({
  obsAnnotations: state.world?.obsAnnotations,
  colorAccessor: state.colors.colorAccessor,
  colorScale: state.colors.scale,
  schema: state.world?.schema,
}))
class Continuous extends React.PureComponent {
  static renderIsStillLoading(zebra, key) {
    return (
      <div
        key={key}
        style={{
          padding: globals.leftSidebarSectionPadding,
          backgroundColor: zebra % 2 === 0 ? globals.lightestGrey : "white",
        }}
      >
        <div
          style={{
            display: "flex",
            justifyContent: "space-between",
            justifyItems: "center",
            alignItems: "center",
          }}
        >
          <div style={{ minWidth: 30 }} />
          <div style={{ display: "flex", alignSelf: "center" }}>
            <span style={{ fontStyle: "italic" }}>{key}</span>
          </div>
          <div
            style={{
              display: "flex",
              justifyContent: "flex-end",
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
      .filter((col) => col.type === "int32" || col.type === "float32")
      .filter((col) => col.name !== obsIndex)
      .map((col) => col.name);

    /* initial value for iterator to simulate index, ranges is an object */
    let zebra = 0;

    return (
      <div>
        {allContinuousNames.map((key) => {
          if (!obsAnnotations.hasCol(key)) {
            // still loading!
            zebra += 1;
            return Continuous.renderIsStillLoading(zebra, key);
          }

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
