// jshint esversion: 6
import _ from "lodash";
import React from "react";
import { connect } from "react-redux";
import Categorical from "./categorical/categorical";
import Continuous from "./continuous/continuous";
import GeneExpression from "./geneExpression";
import ExpressionButtons from "./expression/expressionButtons";
import * as globals from "../globals";
import DynamicScatterplot from "./scatterplot/scatterplot";

@connect(state => ({
  responsive: state.responsive,
  datasetTitle: _.get(state.config, "displayNames.dataset"),
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor
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
    const {
      responsive,
      datasetTitle,
      scatterplotXXaccessor,
      scatterplotYYaccessor
    } = this.props;

    /*
    this magic number should be made less fragile,
    if cellxgene logo or tabs change, this must as well
    */
    const metadataSectionPadding =
      scatterplotXXaccessor && scatterplotYYaccessor ? 450 : 0;
    return (
      <div style={{ position: "fixed" }}>
        <p
          style={{
            position: "fixed",
            left: responsive.width / 2,
            top: 14,
            margin: 0,
            fontSize: 24,
            color: globals.darkerGrey,
            fontWeight: 700,
            width: "100%"
          }}
        >
          cellxgene: &nbsp;
          {datasetTitle}
        </p>

        <div
          style={{
            height: responsive.height - metadataSectionPadding,
            width: 400,
            padding: "0px 10px 10px 10px",
            overflowY: "auto",
            overflowX: "hidden"
          }}
        >
          <Categorical />
          <GeneExpression />
          <Continuous />
        </div>
        <div
          style={{
            position: "fixed",
            bottom: 0,
            left: 0
          }}
        >
          {scatterplotXXaccessor && scatterplotYYaccessor ? (
            <DynamicScatterplot />
          ) : null}
        </div>
        <div style={{ position: "fixed", bottom: 0, right: 0 }}>
          {currentTab === "metadata" ? <ExpressionButtons /> : null}
        </div>
      </div>
    );
  }
}

export default LeftSideBar;
