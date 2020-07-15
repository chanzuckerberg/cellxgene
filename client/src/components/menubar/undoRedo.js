import React from "react";
import { AnchorButton, Tooltip } from "@blueprintjs/core";
import { tooltipHoverOpenDelay } from "../../globals";
import styles from "./menubar.css";

const UndoRedo = React.memo((props) => {
  const { undoDisabled, redoDisabled, dispatch } = props;
  return (
    <div className={`bp3-button-group ${styles.menubarButton}`}>
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
            cursor: "pointer",
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
            cursor: "pointer",
          }}
          data-testid="redo"
        />
      </Tooltip>
    </div>
  );
});

export default UndoRedo;
