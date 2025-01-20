import React from "react";
import { connect } from "react-redux";

import * as globals from "../../globals";
import Logo from "../framework/logo";
import Truncate from "../util/truncate";
import InformationMenu from "./infoMenu";

const DATASET_TITLE_FONT_SIZE = 14;

@connect((state) => {
  const { corpora_props: corporaProps } = state.config;
  const correctVersion =
    ["1.0.0", "1.1.0"].indexOf(corporaProps?.version?.corpora_schema_version) >
    -1;
  return {
    datasetTitle: state.config?.displayNames?.dataset ?? "",
    libraryVersions: state.config?.library_versions,
    aboutLink: state.config?.links?.["about-dataset"],
    tosURL: state.config?.parameters?.about_legal_tos,
    privacyURL: state.config?.parameters?.about_legal_privacy,
    title: correctVersion ? corporaProps?.title : undefined,
  };
})
class LeftSideBar extends React.Component {
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
            Dandelion
          </span>
        </div>
        <div style={{ marginRight: 5, height: "100%" }}>
          <span
            minimal
            style={{
              fontSize: DATASET_TITLE_FONT_SIZE,
              padding: "5px 10px",
            }}
          >
            <Truncate>
              <span style={{ maxWidth: 155 }} data-testid="header">
                {title ?? datasetTitle}
              </span>
            </Truncate>
          </span>
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
