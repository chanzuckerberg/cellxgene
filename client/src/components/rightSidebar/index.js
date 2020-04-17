// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import Continuous from "../continuous/continuous";
import GeneExpression from "../geneExpression";
import * as globals from "../../globals";

@connect(state => ({
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
}))
class RightSidebar extends React.Component {

  render() {
    const { style } = this.props;
    return (
      <div
        ref={ref => { this.ref = ref; }}
        style={{
          /* x y blur spread color */
          borderLeft: `1px solid ${globals.lightGrey}`,
          display: "flex",
          flexDirection: "column",
          position: "relative",
          right: 0,
          top: 0,
          overflowY: "auto",
          ...style
        }}
      >
        <GeneExpression />
        <Continuous />
      </div>
    );
  }
}

export default RightSidebar;
