import React from "react";
import { connect } from "react-redux";
import Categorical from "../categorical";
import * as globals from "../../globals";
import DynamicScatterplot from "../scatterplot/scatterplot";
import TopLeftLogoAndTitle from "./topLeftLogoAndTitle";
import Continuous from "../continuous/continuous";
import GeneExpression from "../geneExpression";

@connect((state) => ({
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
}))
class LeftSideBar extends React.Component {
  render() {
    const { scatterplotXXaccessor, scatterplotYYaccessor } = this.props;
    return (
      <div
        style={{
          /* x y blur spread color */
          borderRight: `1px solid ${globals.lightGrey}`,
          display: "flex",
          flexDirection: "column",
          height: "100%",
        }}
      >
        <TopLeftLogoAndTitle />
        <div
          style={{
            height: "100%",
            width: globals.leftSidebarWidth,
            overflowY: "auto",

            display: "flex",
            flexDirection: "column",
            alignItems: "stretch",
            rowGap: "20px",
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
