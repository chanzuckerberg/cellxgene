// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import Categorical from "./categorical/categorical";
import Continuous from "./continuous/continuous";
import GeneExpression from "./geneExpression";
import * as globals from "../globals";
import DynamicScatterplot from "./scatterplot/scatterplot";
import Logo from "./framework/logo.js";

@connect(state => ({
  responsive: state.responsive,
  datasetTitle: state.config?.displayNames?.dataset,
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
          style={{
            position: "fixed",
            top: globals.cellxgeneTitleTopPadding,
            left: globals.leftSidebarWidth + globals.cellxgeneTitleLeftPadding,
            margin: 0
          }}
        >
          <Logo size={32} />
          <span
            style={{
              fontSize: 28,
              position: "relative",
              top: -4,
              fontWeight: "bold",
              marginLeft: 5,
              color: globals.logoColor
            }}
          >
            cellxgene
          </span>
          <span
            data-testid="header"
            style={{
              fontSize: 16,
              display: "block",
              position: "relative",
              marginTop: 10,
              top: -4
            }}
          >
            {datasetTitle}
          </span>
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
