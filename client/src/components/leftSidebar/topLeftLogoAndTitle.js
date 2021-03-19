import React from "react";
import { connect } from "react-redux";
import { Button } from "@blueprintjs/core";

import * as globals from "../../globals";
import Logo from "../framework/logo";
import Truncate from "../util/truncate";
import InfoDrawer from "../infoDrawer/infoDrawer";
import InformationMenu from "./infoMenu";

const DATASET_TITLE_FONT_SIZE = 14;

@connect((state) => {
  const { corpora_props: corporaProps } = state.config;
  const correctVersion =
    ["1.0.0", "1.1.0"].indexOf(
      corporaProps?.version?.["corpora_schema_version"]
    ) > -1;
  return {
    datasetTitle: state.config?.displayNames?.dataset ?? "",
    libraryVersions: state.config?.["library_versions"],
    aboutLink: state.config?.links?.["about-dataset"],
    tosURL: state.config?.parameters?.["about_legal_tos"],
    privacyURL: state.config?.parameters?.["about_legal_privacy"],
    title: correctVersion ? corporaProps?.title : undefined,
  };
})
class LeftSideBar extends React.Component {
  handleClick = () => {
    const { dispatch } = this.props;
    dispatch({ type: "toggle dataset drawer" });
  };

  render() {
    const {
      datasetTitle,
      libraryVersions,
      aboutLink,
      privacyURL,
      tosURL,
      dispatch,
      title,
    } = this.props;

    return (
      <div
        style={{
          paddingLeft: 8,
          paddingTop: 8,
          width: globals.leftSidebarWidth,
          zIndex: 1,
          borderBottom: `1px solid ${globals.lighterGrey}`,
          display: "flex",
          justifyContent: "space-between",
          alignItems: "center",
        }}
      >
        <div>
          <Logo size={28} />
          <span
            style={{
              fontSize: 24,
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
        </div>
        <div style={{ marginRight: 5, height: "100%" }}>
          <Button
            minimal
            style={{
              fontSize: DATASET_TITLE_FONT_SIZE,
              position: "relative",
              top: -1,
            }}
            onClick={this.handleClick}
          >
            <Truncate>
              <span style={{ maxWidth: 155 }} data-testid="header">
                {title ?? datasetTitle}
              </span>
            </Truncate>
          </Button>
          <InfoDrawer />
          <InformationMenu
            {...{
              libraryVersions,
              aboutLink,
              tosURL,
              privacyURL,
              dispatch,
            }}
          />
        </div>
      </div>
    );
  }
}

export default LeftSideBar;
