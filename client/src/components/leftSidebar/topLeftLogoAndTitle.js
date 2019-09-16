// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import * as globals from "../../globals";
import Logo from "../framework/logo";

@connect(state => ({
  responsive: state.responsive,
  datasetTitle: state.config?.displayNames?.dataset ?? "",
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor
}))
class LeftSideBar extends React.Component {
  render() {
    const { datasetTitle } = this.props;

    const paddingToAvoidScrollBar = 15;

    return (
      <div
        style={{
          paddingLeft: 8,
          paddingTop: 8,
          width: globals.leftSidebarWidth - paddingToAvoidScrollBar,
          position: "absolute",
          backgroundColor: "white",
          zIndex: 8888
          /* x y blur spread color */
          // boxShadow: "-5px -1px 4px 2px rgba(225,225,225,0.4)"
        }}
      >
        <Logo size={30} />
        <span
          style={{
            fontSize: 28,
            position: "relative",
            top: -6,
            fontWeight: "bold",
            marginLeft: 5,
            color: globals.logoColor,
            userSelect: "none"
          }}
        >
          cell
          <span
            style={{
              position: "relative",
              top: 1,
              fontWeight: 300,
              fontSize: 24
            }}
          >
            ×
          </span>
          gene
        </span>
        <div
          data-testid="header"
          style={{
            fontSize: 14,
            position: "relative",
            top: -2,
            display: "inline-block",
            width: "190px",
            marginLeft: "7px",
            height: "1.1em",
            overflow: "hidden",
            wordBreak: "break-all"
          }}
          title={datasetTitle}
        >
          {datasetTitle.length > 29
            ? `${datasetTitle.substring(0, 14)}…${datasetTitle.slice(-14)}`
            : datasetTitle}
        </div>
      </div>
    );
  }
}

export default LeftSideBar;
