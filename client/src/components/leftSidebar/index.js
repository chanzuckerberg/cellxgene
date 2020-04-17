// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import Categorical from "../categorical";
import * as globals from "../../globals";
import DynamicScatterplot from "../scatterplot/scatterplot";
import TopLeftLogoAndTitle from "./topLeftLogoAndTitle";

@connect((state) => ({
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
}))
class LeftSideBar extends React.Component {

  render() {
    const { scatterplotXXaccessor, scatterplotYYaccessor, style } = this.props;
    return (
      <div
        style={{
          /* x y blur spread color */
          borderRight: `1px solid ${globals.lightGrey}`,
          display: "flex",
          flexDirection: "column",
          overflowY: "auto",
          ...style,
        }}
      >
        <TopLeftLogoAndTitle />
        <div
          style={{
            height: "100%",
            width: globals.leftSidebarWidth,
            overflowY: "auto",
          }}
        >
          <Categorical />
        </div>
        {scatterplotXXaccessor && scatterplotYYaccessor ? (
          <DynamicScatterplot />
        ) : null}
      </div>
    );
  }
}

export default LeftSideBar;
