// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import Continuous from "../continuous/continuous";
import GeneExpression from "../geneExpression";
import * as globals from "../../globals";

@connect((state) => ({
  responsive: state.responsive,
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
}))
class RightSidebar extends React.Component {
  render() {
    const { responsive } = this.props;

    return (
      <div
        style={{
          position: "fixed",
          right: 0,
          backgroundColor: "white",
          /* x y blur spread color */
          borderLeft: `1px solid ${globals.lightGrey}`,
        }}
      >
        <div
          style={{
            height: responsive.height,
            width: globals.leftSidebarWidth,
            overflowY: "auto",
            overflowX: "hidden",
          }}
        >
          <GeneExpression />
          <Continuous />
        </div>
      </div>
    );
  }
}

export default RightSidebar;
