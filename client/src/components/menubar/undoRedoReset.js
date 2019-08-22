// jshint esversion: 6
import React from "react";
import { AnchorButton, Tooltip } from "@blueprintjs/core";

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
    <div style={{ marginLeft: 10 }} className="bp3-button-group">
      <Tooltip content="Undo" position="bottom">
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
        />
      </Tooltip>
      <Tooltip content="Redo" position="bottom">
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
        />
      </Tooltip>
      <Tooltip
        content="Reset cellxgene, clearing all selections"
        position="bottom"
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
        >
          reset
        </AnchorButton>
      </Tooltip>
    </div>
  );
}

export default InformationMenu;
