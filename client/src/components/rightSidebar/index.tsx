import React from "react";
import { connect } from "react-redux";
import GeneExpression from "../geneExpression";
import * as globals from "../../globals";

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  scatterplotXXaccessor: (state as any).controls.scatterplotXXaccessor,
  scatterplotYYaccessor: (state as any).controls.scatterplotYYaccessor,
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
          padding: globals.leftSidebarSectionPadding,
        }}
      >
        <GeneExpression />
      </div>
    );
  }
}

export default RightSidebar;
