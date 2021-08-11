// Core dependencies
import { SKELETON } from "@blueprintjs/core/lib/esnext/common/classes";
import React from "react";

// App dependencies
import LeftSidebarSkeleton from "../leftSidebar/leftSidebarSkeleton";
import Layout from "./layout";
import RightSidebarSkeleton from "../rightSidebar/rightSidebarSkeleton";

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
function Skeleton() {
  /*
    Skeleton layout component displayed when in loading state.
    TODO(cc)
     - Remove dupe of "graph" area inline styles
    */
  return (
    <Layout>
      <LeftSidebarSkeleton />
      {() => (
        <>
          <div
            style={{
              display: "flex",
              justifyContent: "space-between",
              left: 8,
              position: "absolute",
              right: 8,
              top: 8,
              zIndex: 3,
            }}
          >
            <div
              style={{ height: "30px", width: "calc(100% - 482px - 10px)" }}
              className={SKELETON}
            />
            <div
              style={{ height: "30px", width: "482px" }}
              className={SKELETON}
            />
          </div>
        </>
      )}
      <RightSidebarSkeleton />
    </Layout>
  );
}

export default Skeleton;
