// jshint esversion: 6
import React from "react";
import _ from "lodash";
import Categorical from "./categorical/categorical";
import Continuous from "./continuous/continuous";
import ExpressionButtons from "./expression/expressionButtons";
import { connect } from "react-redux";
import Heatmap from "./expression/diffExpHeatmap";
import * as globals from "../globals";
import DynamicScatterplot from "./scatterplot/scatterplot";

@connect(state => {
  return {
    responsive: state.responsive
  };
})
class LeftSideBar extends React.Component {
  constructor(props) {
    super(props);
    this.metadataSectionPadding = 300;
    this.state = {
      currentTab: "metadata"
    };
  }
  render() {
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
          CELLxGENE {globals.datasetTitle}{" "}
        </p>
        <div style={{ padding: 10 }}>
          <button
            style={{
              padding: "none",
              outline: 0,
              fontSize: 14,
              fontWeight: this.state.currentTab === "metadata" ? 700 : 400,
              fontStyle:
                this.state.currentTab === "metadata" ? "inherit" : "italic",
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
            style={{
              padding: "none",
              outline: 0,
              fontSize: 14,
              fontWeight: this.state.currentTab === "expression" ? 700 : 400,
              fontStyle:
                this.state.currentTab === "expression" ? "inherit" : "italic",
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
            height: this.props.responsive.height - this.metadataSectionPadding,
            width: 400,
            padding: 10,
            overflowY: "auto",
            overflowX: "hidden"
          }}
        >
          {this.state.currentTab === "metadata" ? <Categorical /> : null}
          {this.state.currentTab === "metadata" ? <Continuous /> : null}
          {this.state.currentTab === "expression" ? <Heatmap /> : null}
        </div>
        <div style={{ position: "fixed", bottom: 0, left: 0 }}>
          {this.state.currentTab === "metadata" ? <ExpressionButtons /> : null}
          {this.state.currentTab === "expression" ? (
            <DynamicScatterplot />
          ) : null}
        </div>
      </div>
    );
  }
}

export default LeftSideBar;
