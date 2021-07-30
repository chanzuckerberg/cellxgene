import React from "react";
import { AnchorButton, ButtonGroup, Tooltip } from "@blueprintjs/core";
// @ts-expect-error ts-migrate(2307) FIXME: Cannot find module './menubar.css' or its correspo... Remove this comment to see the full error message
import styles from "./menubar.css";
import * as globals from "../../globals";

const Subset = React.memo((props) => {
  const {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'subsetPossible' does not exist on type '... Remove this comment to see the full error message
    subsetPossible,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'subsetResetPossible' does not exist on t... Remove this comment to see the full error message
    subsetResetPossible,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'handleSubset' does not exist on type '{ ... Remove this comment to see the full error message
    handleSubset,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'handleSubsetReset' does not exist on typ... Remove this comment to see the full error message
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
