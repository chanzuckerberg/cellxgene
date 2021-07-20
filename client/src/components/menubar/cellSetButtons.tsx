import React from "react";
import { AnchorButton, Tooltip } from "@blueprintjs/core";
import { connect } from "react-redux";
import { tooltipHoverOpenDelay } from "../../globals";
import actions from "../../actions";

@connect((state) => ({
  differential: state.differential,
}))
class CellSetButton extends React.PureComponent {
  set() {
    const { dispatch, eitherCellSetOneOrTwo } = this.props;

    dispatch(actions.setCellSetFromSelection(eitherCellSetOneOrTwo));
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
