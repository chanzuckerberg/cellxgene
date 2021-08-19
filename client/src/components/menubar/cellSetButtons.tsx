import React from "react";
import { AnchorButton, Tooltip } from "@blueprintjs/core";
import { connect } from "react-redux";
import { tooltipHoverOpenDelay } from "../../globals";
import actions from "../../actions";

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  differential: (state as any).differential,
}))
class CellSetButton extends React.PureComponent {
  set() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, eitherCellSetOneOrTwo } = this.props;
    dispatch(actions.setCellSetFromSelection(eitherCellSetOneOrTwo));
  }

  render() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'differential' does not exist on type 'Re... Remove this comment to see the full error message
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
