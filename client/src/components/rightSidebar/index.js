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
        <Chat />
      </div>
    );
  }
}

export default RightSidebar;
