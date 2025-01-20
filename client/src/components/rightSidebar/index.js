import React from "react";
import { connect } from "react-redux";
import * as globals from "../../globals";

@connect((state) => ({
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
}))
class RightSidebar extends React.Component {
  render() {
    return (
      <div
        style={{
          /* x y blur spread color */
          borderLeft: `1px solid ${globals.lightGrey}`,
          display: "flex",
          flexDirection: "column",
          position: "relative",
          overflowY: "inherit",
          height: "inherit",
          width: "inherit",
        }}
      >
        <div
          style={{
            flex: 1,
            borderBottom: `0.5px solid ${globals.lightGrey}`,
            overflowY: "inherit",
            padding: globals.leftSidebarSectionPadding,
          }}
        />
        <div
          style={{
            flex: 1,
            borderTop: `0.5px solid ${globals.lightGrey}`,
            overflowY: "inherit",
            padding: globals.leftSidebarSectionPadding,
          }}
        />
      </div>
    );
  }
}

export default RightSidebar;
