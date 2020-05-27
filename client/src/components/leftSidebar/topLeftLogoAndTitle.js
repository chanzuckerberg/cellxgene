// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import * as globals from "../../globals";
import Logo from "../framework/logo";
import Truncate from "../util/truncate";

@connect((state) => ({
  datasetTitle: state.config?.displayNames?.dataset ?? "",
  aboutURL: state.config?.links?.["about-dataset"],
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
}))
class LeftSideBar extends React.Component {
  render() {
    const { datasetTitle, aboutURL } = this.props;
    const width = 190;
    const fontSize = 14;

    return (
      <div
        style={{
          paddingLeft: 8,
          paddingTop: 8,
          width: globals.leftSidebarWidth,
          zIndex: 1,
          borderBottom: `1px solid ${globals.lighterGrey}`,
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
            userSelect: "none",
          }}
        >
          cell
          <span
            style={{
              position: "relative",
              top: 1,
              fontWeight: 300,
              fontSize: 24,
            }}
          >
            ×
          </span>
          gene
        </span>
        <div
          data-testid="header"
          style={{
            fontSize,
            position: "relative",
            top: -6,
            display: "inline-block",
            width,
            marginLeft: "7px",
            height: "1.2em",
            overflow: "hidden",
            wordBreak: "break-all",
          }}
        >
          {aboutURL ? (
            <Truncate size={width} fontSize={fontSize}>
              <a href={aboutURL} target="_blank" rel="noopener noreferrer">
                {datasetTitle}
              </a>
            </Truncate>
          ) : (
            <Truncate size={width} fontSize={fontSize}>
              <span>{datasetTitle}</span>
            </Truncate>
          )}
        </div>
      </div>
    );
  }
}

export default LeftSideBar;
