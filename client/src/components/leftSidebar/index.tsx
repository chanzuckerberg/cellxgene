import React from "react";
import { connect } from "react-redux";
import Categorical from "../categorical";
import * as globals from "../../globals";
import DynamicScatterplot from "../scatterplot/scatterplot";
import TopLeftLogoAndTitle from "./topLeftLogoAndTitle";
import Continuous from "../continuous/continuous";

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  scatterplotXXaccessor: (state as any).controls.scatterplotXXaccessor,
  scatterplotYYaccessor: (state as any).controls.scatterplotYYaccessor,
}))
class LeftSideBar extends React.Component {
  render() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'scatterplotXXaccessor' does not exist on... Remove this comment to see the full error message
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
          }}
        >
          <Categorical />
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
