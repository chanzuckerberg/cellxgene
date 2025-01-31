import React from "react";
import * as globals from "../../globals";

import Chat from "../chat";

class RightSidebar extends React.PureComponent {
  render() {
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
        }}
      >
        <div
          style={{
            flex: 1,
            borderBottom: `0.5px solid ${globals.lightGrey}`,
            overflowY: "inherit",
            padding: globals.leftSidebarSectionPadding,
          }}
        />
        <Chat />
      </div>
    );
  }
}

export default RightSidebar;
