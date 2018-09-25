// jshint esversion: 6
import _ from "lodash";
import React from "react";
import { connect } from "react-redux";
import Categorical from "./categorical/categorical";
import Continuous from "./continuous/continuous";
import ExpressionButtons from "./expression/expressionButtons";
import Heatmap from "./expression/diffExpHeatmap";
import * as globals from "../globals";
import DynamicScatterplot from "./scatterplot/scatterplot";

@connect(state => ({
  responsive: state.responsive,
  datasetTitle: _.get(state.config, "displayNames.dataset")
}))
class LeftSideBar extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      currentTab: "metadata"
    };
  }

  render() {
    const { currentTab } = this.state;
    const { responsive, datasetTitle } = this.props;

    /*
    this magic number should be made less fragile,
    if cellxgene logo or tabs change, this must as well
    */
    const metadataSectionPadding = currentTab === "metadata" ? 88 : 450;
    return (
      <div style={{ position: "fixed" }}>
        <p
          style={{
            margin: 10,
            fontSize: 24,
            color: globals.lightGrey,
            fontWeight: 700,
            width: "100%"
          }}
        >
          cellxgene &nbsp;
          {datasetTitle}
        </p>
        <div style={{ padding: 10 }}>
          <button
            type="button"
            style={{
              padding: "none",
              outline: 0,
              fontSize: 14,
              fontWeight: currentTab === "metadata" ? 700 : 400,
              fontStyle: currentTab === "metadata" ? "inherit" : "italic",
              cursor: "pointer",
              border: "none",
              backgroundColor: "#FFF",
              borderTop: "none",
              borderBottom: "none",
              borderRight: "none",
              borderLeft: "none"
            }}
            onClick={() => {
              this.setState({ currentTab: "metadata" });
            }}
          >
            Metadata
          </button>
          <button
            type="button"
            style={{
              padding: "none",
              outline: 0,
              fontSize: 14,
              fontWeight: currentTab === "expression" ? 700 : 400,
              fontStyle: currentTab === "expression" ? "inherit" : "italic",
              cursor: "pointer",
              border: "none",
              backgroundColor: "#FFF",
              borderTop: "none",
              borderBottom: "none",
              borderRight: "none",
              borderLeft: "none"
            }}
            onClick={() => {
              this.setState({ currentTab: "expression" });
            }}
          >
            Expression
          </button>
        </div>
        <div
          style={{
            height: responsive.height - metadataSectionPadding,
            width: 400,
            padding: 10,
            overflowY: "auto",
            overflowX: "hidden"
          }}
        >
          {currentTab === "metadata" ? <Categorical /> : null}
          {currentTab === "metadata" ? <Continuous /> : null}
          {currentTab === "expression" ? <Heatmap /> : null}
        </div>
        <div
          style={{
            position: "fixed",
            bottom: 0,
            left: 0
          }}
        >
          {currentTab === "expression" ? <DynamicScatterplot /> : null}
        </div>
        <div style={{ position: "fixed", bottom: 0, right: 0 }}>
          {currentTab === "metadata" ? <ExpressionButtons /> : null}
        </div>
      </div>
    );
  }
}

export default LeftSideBar;
