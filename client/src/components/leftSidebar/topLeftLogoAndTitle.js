// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import * as globals from "../../globals";
import Logo from "../framework/logo";
import Truncate from "../util/truncate";
import InfoDrawer from "../infoDrawer/infoDrawer";

const DATASET_TITLE_WIDTH = 190;
const DATASET_TITLE_FONT_SIZE = 14;

@connect((state) => ({
  datasetTitle: state.config?.displayNames?.dataset ?? "",
  aboutURL: state.config?.links?.["about-dataset"],
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
}))
class LeftSideBar extends React.Component {
  render() {
    const { datasetTitle, aboutURL } = this.props;

    return (
      <InfoDrawer title="Dataset Info" {...{ aboutURL, datasetTitle }}>
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
            style={{
              fontSize: DATASET_TITLE_FONT_SIZE,
              position: "relative",
              top: -6,
              display: "inline-block",
              width: DATASET_TITLE_WIDTH,
              marginLeft: "7px",
              height: "1.2em",
              overflow: "hidden",
              wordBreak: "break-all",
            }}
          >
            {aboutURL ? (
              <Truncate>
                <span style={{ width: 185 }} data-testid="header">
                  {datasetTitle}
                </span>
              </Truncate>
            ) : (
              <Truncate>
                <span style={{ width: 185 }} data-testid="header">
                  {datasetTitle}
                </span>
              </Truncate>
            )}
          </div>
        </div>
      </InfoDrawer>
    );
  }
}

export default LeftSideBar;
