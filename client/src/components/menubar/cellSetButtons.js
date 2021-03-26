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
    const { differential, dispatch, eitherCellSetOneOrTwo } = this.props;

    if (!differential.diffExp) {
      // disallow this action if the user has active differential expression results
      dispatch(actions.setCellSetFromSelection(eitherCellSetOneOrTwo));
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
