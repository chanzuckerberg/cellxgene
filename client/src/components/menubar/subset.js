import React from "react";
import {AnchorButton, ButtonGroup, Tooltip} from "@blueprintjs/core";
import * as globals from "../../globals";


function Subset(props) {
  const {
    subsetPossible,
    subsetResetPossible,
    handleSubset,
    handleSubsetReset,
  } = props;

  return (
    <ButtonGroup style={{marginRight: "10px"}}>
      <Tooltip
        content="Subset to currently selected cells and associated metadata"
        position="bottom"
        hoverOpenDelay={globals.tooltipHoverOpenDelay}
      >
        <AnchorButton
          data-testid="subset-button"
          active={!subsetPossible}
          icon="pie-chart"
          onClick={() => {
            if (subsetPossible) handleSubset();
          }}
        />
      </Tooltip>
      <Tooltip
        content="Undo subset and show all cells and associated metadata"
        position="bottom"
        hoverOpenDelay={globals.tooltipHoverOpenDelay}
      >
        <AnchorButton
          data-testid="reset-subset-button"
          active={!subsetResetPossible}
          icon="full-circle"
          onClick={() => {
            if (subsetResetPossible) handleSubsetReset();
          }}
        />
      </Tooltip>
    </ButtonGroup>
  );
}

export default Subset;
