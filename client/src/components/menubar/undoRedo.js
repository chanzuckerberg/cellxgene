import React from "react";
import { AnchorButton, ButtonGroup, Tooltip } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { tooltipHoverOpenDelay } from "../../globals";
import styles from "./menubar.css";

const UndoRedo = React.memo((props) => {
  const { undoDisabled, redoDisabled, dispatch } = props;
  return (
    <ButtonGroup className={`${styles.menubarButton}`}>
      <Tooltip
        content="Undo"
        position="bottom"
        hoverOpenDelay={tooltipHoverOpenDelay}
      >
        <AnchorButton
          type="button"
          icon={IconNames.UNDO}
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
          icon={IconNames.REDO}
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
    </ButtonGroup>
  );
});

export default UndoRedo;
