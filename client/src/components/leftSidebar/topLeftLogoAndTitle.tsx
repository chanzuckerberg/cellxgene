import React from "react";
import { connect } from "react-redux";

import * as globals from "../../globals";
import Logo from "../framework/logo";
import Title from "../framework/title";
import InformationMenu from "./infoMenu";

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  libraryVersions: (state as any).config?.library_versions,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  aboutLink: (state as any).config?.links?.["about-dataset"],
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  tosURL: (state as any).config?.parameters?.about_legal_tos,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  privacyURL: (state as any).config?.parameters?.about_legal_privacy,
}))
class LeftSideBar extends React.Component {
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'libraryVersions' does not exist on type ... Remove this comment to see the full error message
      libraryVersions,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'aboutLink' does not exist on type 'Reado... Remove this comment to see the full error message
      aboutLink,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'privacyURL' does not exist on type 'Read... Remove this comment to see the full error message
      privacyURL,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'tosURL' does not exist on type 'Readonly... Remove this comment to see the full error message
      tosURL,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
      dispatch,
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
          <Title />
        </div>
        <div style={{ marginRight: 5, height: "100%" }}>
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
