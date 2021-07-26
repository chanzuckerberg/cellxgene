import React from "react";
import { AnchorButton, ButtonGroup, Tooltip } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { tooltipHoverOpenDelay } from "../../globals";
// @ts-expect-error ts-migrate(2307) FIXME: Cannot find module './menubar.css' or its correspo... Remove this comment to see the full error message
import styles from "./menubar.css";

const UndoRedo = React.memo((props) => {
  // @ts-expect-error ts-migrate(2339) FIXME: Property 'undoDisabled' does not exist on type '{ ... Remove this comment to see the full error message
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
