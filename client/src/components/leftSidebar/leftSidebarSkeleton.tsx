// Core dependencies
import { SKELETON } from "@blueprintjs/core/lib/esnext/common/classes";
import React from "react";

// App dependencies
import Logo from "../framework/logo";
import Title from "../framework/title";
import * as globals from "../../globals";

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
function LeftSidebarSkeleton() {
  /*
    Skeleton of left side bar, to be displayed during data load.
    TODO(cc)
      - Remove dupe of LeftSidebar inline styles
      - Remove dupe of TopLeftLogoAndTitle styles
     */
  return (
    <div
      style={{
        borderRight: `1px solid ${globals.lightGrey}`,
        display: "flex",
        flexDirection: "column",
        height: "100%",
      }}
    >
      {/* TopLeftLogoAndTitle */}
      <div
        style={{
          paddingLeft: 8,
          paddingRight: 5,
          paddingTop: 8,
          width: globals.leftSidebarWidth,
          zIndex: 1,
          borderBottom: `1px solid ${globals.lighterGrey}`,
          display: "flex",
          justifyContent: "space-between",
          alignItems: "flex-start",
        }}
      >
        <div>
          <Logo size={28} />
          <Title />
        </div>
        {/* Hamburger */}
        <div style={{ height: 30, width: 30 }} className={SKELETON} />
      </div>
      {/* Categorical */}
      <div style={{ padding: 8 }}>
        {[...Array(10).keys()].map((i) => (
          <div
            key={i}
            style={{ height: 30, marginBottom: 4 }}
            className={SKELETON}
          />
        ))}
      </div>
      {/* Continuous */}
      {[...Array(2).keys()].map((i) => (
        <div key={i} style={{ height: 211 }} className={SKELETON} />
      ))}
    </div>
  );
}

export default LeftSidebarSkeleton;
