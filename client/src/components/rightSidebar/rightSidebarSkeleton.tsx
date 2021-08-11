// Core dependencies
import { SKELETON } from "@blueprintjs/core/lib/esnext/common/classes";
import React from "react";

// App dependencies
import * as globals from "../../globals";

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
function RightSidebarSkeleton() {
  /*
  Skeleton of left side bar, to be displayed during data load.
  TODO(cc)
    - Remove dupe of RightSidebar inline styles
 */
  return (
    <div
      style={{
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
      {/* Create new gene set button */}
      <div
        style={{
          marginBottom: 10,
          position: "relative",
          top: -2,
          height: "30px",
          width: "133px",
        }}
        className={SKELETON}
      />
    </div>
  );
}

export default RightSidebarSkeleton;
