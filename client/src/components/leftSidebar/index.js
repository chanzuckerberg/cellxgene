// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import Categorical from "../categorical/categorical";
import Continuous from "../continuous/continuous";
import GeneExpression from "../geneExpression";
import * as globals from "../../globals";
import DynamicScatterplot from "../scatterplot/scatterplot";
import TopLeftLogoAndTitle from "./topLeftLogoAndTitle";

@connect(state => ({
  responsive: state.responsive,
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor
}))
class LeftSideBar extends React.Component {
  render() {
    const {
      responsive,
      scatterplotXXaccessor,
      scatterplotYYaccessor
    } = this.props;

    /*
    this magic number should be made less fragile,
    if cellxgene logo or tabs change, this must as well
    */
    const metadataSectionPadding = 0;
    const logoRelatedPadding = 50;

    return (
      <div
        style={{
          position: "fixed",
          backgroundColor: "white",
          /* x y blur spread color */
          boxShadow: "-3px 0px 6px 2px rgba(153,153,153,0.4)"
        }}
      >
        <TopLeftLogoAndTitle />

        <div
          style={{
            height: responsive.height - metadataSectionPadding,
            paddingTop: logoRelatedPadding,
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
