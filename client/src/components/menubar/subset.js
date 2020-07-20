import React from "react";
import { AnchorButton, ButtonGroup, Tooltip } from "@blueprintjs/core";
import styles from "./menubar.css";
import * as globals from "../../globals";

const Subset = React.memo((props) => {
  const {
    subsetPossible,
    subsetResetPossible,
    handleSubset,
    handleSubsetReset,
  } = props;

  return (
    <ButtonGroup className={styles.menubarButton}>
      <Tooltip
        content="Subset to currently selected cells and associated metadata"
        position="bottom"
        hoverOpenDelay={globals.tooltipHoverOpenDelay}
      >
        <AnchorButton
          type="button"
          data-testid="subset-button"
          disabled={!subsetPossible}
          icon="pie-chart"
          onClick={handleSubset}
        />
      </Tooltip>
      <Tooltip
        content="Undo subset and show all cells and associated metadata"
        position="bottom"
        hoverOpenDelay={globals.tooltipHoverOpenDelay}
      >
        <AnchorButton
          type="button"
          data-testid="reset-subset-button"
          disabled={!subsetResetPossible}
          icon="full-circle"
          onClick={handleSubsetReset}
        />
      </Tooltip>
    </ButtonGroup>
  );
});

export default Subset;
