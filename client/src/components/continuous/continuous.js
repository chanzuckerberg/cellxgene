// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import { connect } from "react-redux";
import { Button } from "@blueprintjs/core";

import * as globals from "../../globals";
import HistogramBrush from "../brushableHistogram";

@connect((state) => ({
  schema: state.annoMatrix?.schema,
  annoMatrix: state.annoMatrix,
}))
class Continuous extends React.PureComponent {
  constructor(props) {
    super(props);

    /*
    allContinuous: map of category name to object.  Object contains:
      * status: loading status for the category
      * summary: result of Dataframe.summary()
    */
    this.state = {
      allContinuous: new Map(),
    };
  }

  setAllContinuous(colName, status) {
    this.setState({
      allContinuous: new Map(this.state.allContinuous).set(colName, status),
    });
  }

  fetchState(prevProps) {
    const { schema, annoMatrix } = this.props;
    let allContinuousNames;

    if (!schema) return;
    if (schema !== prevProps?.schema) {
      const obsIndex = schema.annotations.obs.index;
      allContinuousNames = schema.annotations.obs.columns
        .filter((col) => col.type === "int32" || col.type === "float32")
        .filter((col) => col.name !== obsIndex)
        .map((col) => col.name);

      // Map preserves order of insertion, allowing ordering of component render.
      this.setState({
        allContinuous: new Map(
          allContinuousNames.map((name) => [name, { status: "pending" }])
        ),
      });
    }

    if (!annoMatrix) return;
    if (annoMatrix !== prevProps?.annoMatrix) {
      allContinuousNames.forEach((colName) =>
        annoMatrix.fetch("obs", colName).then(
          (df) =>
            this.setAllContinuous(colName, {
              status: "success",
              summary: df.col(colName).summarize(),
            }),
          (error) => this.setAllContinuous(colName, { status: "error" })
        )
      );
    }
  }

  componentDidMount() {
    const { annoMatrix } = this.props;
    if (annoMatrix) this.fetchState();
  }

  componentDidUpdate(prevProps) {
    const { annoMatrix } = this.props;
    if (annoMatrix !== prevProps.annoMatrix) {
      this.fetchState();
    }
  }

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
    const { schema, annoMatrix } = this.props;

    /* initial value for iterator to simulate index, ranges is an object */
    let zebra = 0;

    const { allContinuous } = this.state;
    return (
      <div>
        {Array.from(allContinuous, ([key, val]) => {
          const { status, summary } = val;
          if (status === "error") return null; // skip
          if (status === "pending") {
            // sitll loading!
            zebra += 1;
            return Continuous.renderIsStillLoading(zebra, key);
          }

          const nonFiniteExtent =
            summary.min === undefined ||
            summary.max === undefined ||
            Number.isNaN(summary.min) ||
            Number.isNaN(summary.max);
          if (
            status === "success" &&
            !summary.categorical &&
            !nonFiniteExtent
          ) {
            zebra += 1;
            return (
              <HistogramBrush
                key={key}
                field={key}
                isObs
                zebra={zebra % 2 === 0}
                // ranges={summary}
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
