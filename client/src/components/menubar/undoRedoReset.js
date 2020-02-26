// jshint esversion: 6
import React from "react";
import { AnchorButton, Tooltip } from "@blueprintjs/core";
import { tooltipHoverOpenDelay } from "../../globals";

function InformationMenu(props) {
  const {
    resettingInterface,
    undoDisabled,
    redoDisabled,
    resetInterface,
    isResetDisabled,
    dispatch
  } = props;
  return (
    <div style={{ marginRight: 10 }} className="bp3-button-group">
      <Tooltip
        content="Undo"
        position="bottom"
        hoverOpenDelay={tooltipHoverOpenDelay}
      >
        <AnchorButton
          type="button"
          className="bp3-button bp3-icon-undo"
          disabled={undoDisabled}
          onClick={() => {
            dispatch({ type: "@@undoable/undo" });
          }}
          style={{
            cursor: "pointer"
          }}
          data-testid="undo"
        />
      </Tooltip>
      <Tooltip
        content="Redo"
        position="bottom"
        hoverOpenDelay={tooltipHoverOpenDelay}
      >
        <AnchorButton
          type="button"
          className="bp3-button bp3-icon-redo"
          disabled={redoDisabled}
          onClick={() => {
            dispatch({ type: "@@undoable/redo" });
          }}
          style={{
            cursor: "pointer"
          }}
          data-testid="redo"
        />
      </Tooltip>
      <Tooltip
        content="Reset cellxgene, clearing all selections"
        position="bottom"
        hoverOpenDelay={tooltipHoverOpenDelay}
      >
        <AnchorButton
          disabled={isResetDisabled()}
          style={{ marginLeft: 10 }}
          type="button"
          loading={resettingInterface}
          intent="none"
          icon="refresh"
          onClick={resetInterface}
          data-testid="reset"
          data-testclass={`resetting-${resettingInterface}`}
        />
      </Tooltip>
    </div>
  );
}

export default InformationMenu;
