// jshint esversion: 6
import React from "react";
import { AnchorButton, Tooltip } from "@blueprintjs/core";
import { connect } from "react-redux";
import { tooltipHoverOpenDelay } from "../../globals";

@connect()
class CellSetButton extends React.PureComponent {
  set() {
    const {
      differential,
      crossfilter,
      dispatch,
      eitherCellSetOneOrTwo,
    } = this.props;

    let set = crossfilter.allSelectedLabels();
    if (set.length === 0) set = null;
    if (!differential.diffExp) {
      /* diffexp needs to be cleared before we store a new set */
      dispatch({
        type: `store current cell selection as differential set ${eitherCellSetOneOrTwo}`,
        data: set,
      });
    }
  }

  render() {
    const { differential, eitherCellSetOneOrTwo } = this.props;
    const cellListName = `celllist${eitherCellSetOneOrTwo}`;
    const cellsSelected = differential[cellListName]
      ? differential[cellListName].length
      : 0;
    return (
      <Tooltip
        content="Save current selection for differential expression computation"
        position="bottom"
        hoverOpenDelay={tooltipHoverOpenDelay}
      >
        <AnchorButton
          type="button"
          disabled={differential.diffExp}
          onClick={() => {
            this.set();
          }}
          data-testid={`cellset-button-${eitherCellSetOneOrTwo}`}
        >
          {eitherCellSetOneOrTwo}
          {": "}
          <span data-testid={`cellset-count-${eitherCellSetOneOrTwo}`}>
            {cellsSelected}
          </span>
          {" cells"}
        </AnchorButton>
      </Tooltip>
    );
  }
}

export default CellSetButton;
