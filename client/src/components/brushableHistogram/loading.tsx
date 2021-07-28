import React from "react";
import { Button } from "@blueprintjs/core";

import * as globals from "../../globals";

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
const StillLoading = ({ zebra, displayName }: any) => {
  /*
    Render a loading indicator for the field.
    */
  return (
    <div
      data-testclass="gene-loading-spinner"
      style={{
        padding: globals.leftSidebarSectionPadding,
        backgroundColor: zebra ? globals.lightestGrey : "white",
      }}
    >
      <div
        style={{
          display: "flex",
          justifyContent: "space-between",
          justifyItems: "center",
          alignItems: "center",
        }}
      >
        <div style={{ minWidth: 30 }} />
        <div style={{ display: "flex", alignSelf: "center" }}>
          <span style={{ fontStyle: "italic" }}>{displayName}</span>
        </div>
        <div
          style={{
            display: "flex",
            justifyContent: "flex-end",
          }}
        >
          <Button minimal loading intent="primary" />
        </div>
      </div>
    </div>
  );
};

export default StillLoading;
