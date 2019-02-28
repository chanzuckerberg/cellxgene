// jshint esversion: 6
import _ from "lodash";
import React from "react";
import { connect } from "react-redux";
import Categorical from "./categorical/categorical";
import Continuous from "./continuous/continuous";
import GeneExpression from "./geneExpression";
import * as globals from "../globals";
import DynamicScatterplot from "./scatterplot/scatterplot";

@connect(state => ({
  responsive: state.responsive,
  datasetTitle: _.get(state.config, "displayNames.dataset"),
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor
}))
class LeftSideBar extends React.Component {
  render() {
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
    const metadataSectionPadding = 0;
    // scatterplotXXaccessor && scatterplotYYaccessor ? 450 : 0;

    return (
      <div
        style={{
          position: "fixed",
          backgroundColor: "white",
          /* x y blur spread color */
          boxShadow: "1px 0px 6px 2px rgba(153,153,153,0.4)"
        }}
      >
        <p
          data-testid="header"
          style={{
            position: "fixed",
            top: globals.cellxgeneTitleTopPadding,
            left: globals.leftSidebarWidth + globals.cellxgeneTitleLeftPadding,
            margin: 0,
            fontSize: globals.largestFontSize,
            color: globals.darkerGrey,
            width: "100%"
          }}
        >
          cellxgene: {datasetTitle}
        </p>

        <div
          style={{
            height: responsive.height - metadataSectionPadding,
            width: globals.leftSidebarWidth,
            overflowY: "auto",
            overflowX: "hidden"
          }}
        >
          <Categorical />
          <GeneExpression />
          <Continuous />
        </div>
        {scatterplotXXaccessor && scatterplotYYaccessor ? (
          <DynamicScatterplot />
        ) : null}
      </div>
    );
  }
}

export default LeftSideBar;
