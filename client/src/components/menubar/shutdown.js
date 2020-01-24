// jshint esversion: 6
import React from "react";
import { AnchorButton, Tooltip } from "@blueprintjs/core";
import { tooltipHoverOpenDelay } from "../../globals";

function InformationMenu(props) {
  const {
    shutdownServer,
    dispatch,
    handleShutdownServer
  } = props;
  return (
    <div style={{ marginRight: 10 }} className="bp3-button-group">
      <Tooltip
        content="Shudown server"
        position="bottom"
        hoverOpenDelay={tooltipHoverOpenDelay}
      >
        <AnchorButton
          type="button"
          className="bp3-button bp3-icon-power"
          disabled={!shutdownServer}
          onClick={handleShutdownServer}
          style={{
            cursor: "pointer"
          }}
        />
      </Tooltip>
      
    </div>
  );
}

export default InformationMenu;
