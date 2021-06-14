import React from "react";
import { connect } from "react-redux";

import * as globals from "../../globals";
import Logo from "../framework/logo";
import InformationMenu from "./infoMenu";

@connect((state) => {
  const selectedDatasetId = state.collections?.selectedDatasetId;
  const collection = state.collections?.collectionsByDatasetId?.get(
    selectedDatasetId
  );
  return {
    libraryVersions: state.config?.library_versions,
    aboutLink: state.config?.links?.["about-dataset"],
    dataSourceLink: collection?.links.find(
      (link) => link.type === "DATA_SOURCE" // TODO(cc) constants for link types, remove type from returned link prop?
    ),
    tosURL: state.config?.parameters?.about_legal_tos,
    privacyURL: state.config?.parameters?.about_legal_privacy,
    skeleton: state.skeleton.skeleton,
  };
})
class LeftSideBar extends React.Component {
  render() {
    const {
      libraryVersions,
      aboutLink,
      dataSourceLink,
      privacyURL,
      tosURL,
      dispatch,
      skeleton,
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
          <InformationMenu
            {...{
              libraryVersions,
              aboutLink,
              dataSourceLink,
              tosURL,
              privacyURL,
              dispatch,
              skeleton,
            }}
          />
        </div>
      </div>
    );
  }
}

export default LeftSideBar;
