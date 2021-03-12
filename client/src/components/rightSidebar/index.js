import React from "react";
import { connect } from "react-redux";
import GeneExpression from "../geneExpression";
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
          padding: globals.leftSidebarSectionPadding,
        }}
      >
        {/* <div
          style={{
            border: "1px dotted #738694",
            background: "#F5F8FA",
            width: "100%",
            textAlign: "center",
            marginBottom: 5,
            position: "relative",
            top: -2,
            padding: 5,
          }}
        >
          upload geneset csv
        </div> */}
        <GeneExpression />
      </div>
    );
  }
}

export default RightSidebar;
