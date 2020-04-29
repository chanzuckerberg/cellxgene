// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import {
  Popover,
  Button,
  ButtonGroup,
  AnchorButton,
  Tooltip,
  Position,
} from "@blueprintjs/core";
import * as globals from "../../globals";
import actions from "../../actions";
import CellSetButton from "./cellSetButtons";

@connect((state) => ({
  config: state.config,
  crossfilter: state.crossfilter,
  differential: state.differential,
  celllist1: state.differential?.celllist1,
  celllist2: state.differential?.celllist2,
  diffexpMayBeSlow: state.config?.parameters?.diffexpMayBeSlow ?? false,
  diffexpCellcountMax: state.config?.limits?.diffexpCellcountMax,
}))
class DiffexpButtons extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      userDismissedPopover: false,
    };
  }

  computeDiffExp = () => {
    const { dispatch, differential } = this.props;
    if (differential.celllist1 && differential.celllist2) {
      dispatch(
        actions.requestDifferentialExpression(
          differential.celllist1,
          differential.celllist2
        )
      );
    }
  };

  clearDifferentialExpression = () => {
    const { dispatch, differential } = this.props;
    dispatch({
      type: "clear differential expression",
      diffExp: differential.diffExp,
    });
    dispatch({
      type: "clear scatterplot",
    });
  };

  handlePopoverDismiss = () => {
    this.setState({
      userDismissedPopover: true,
    });
  };

  render() {
    /* diffexp-related buttons may be disabled */
    const { differential, diffexpMayBeSlow, diffexpCellcountMax } = this.props;
    const { userDismissedPopover } = this.state;

    const haveBothCellSets =
      !!differential.celllist1 && !!differential.celllist2;

    const haveEitherCellSet =
      !!differential.celllist1 || !!differential.celllist2;

    const slowMsg = diffexpMayBeSlow
      ? " (CAUTION: large dataset - may take longer or fail)"
      : "";
    const tipMessage = `See top 10 differentially expressed genes${slowMsg}`;
    const tipMessageWarn = `The total number of cells for differential expression computation
                            may not exceed ${diffexpCellcountMax}. Try reselecting new cell sets.`;

    const warnMaxSizeExceeded =
      haveEitherCellSet &&
      !!diffexpCellcountMax &&
      (differential.celllist1?.length ?? 0) +
        (differential.celllist2?.length ?? 0) >
        diffexpCellcountMax;

    return (
      <ButtonGroup style={{ marginRight: 10 }}>
        <CellSetButton
          {...this.props} // eslint-disable-line react/jsx-props-no-spreading
          eitherCellSetOneOrTwo={1}
        />
        <CellSetButton
          {...this.props} // eslint-disable-line react/jsx-props-no-spreading
          eitherCellSetOneOrTwo={2}
        />
        {!differential.diffExp ? (
          <Popover
            isOpen={/* warnMaxSizeExceeded && !userDismissedPopover */ false}
            position={Position.BOTTOM}
            target={
              <Tooltip
                content={warnMaxSizeExceeded ? tipMessageWarn : tipMessage}
                position="bottom"
                hoverOpenDelay={globals.tooltipHoverOpenDelayQuick}
                intent={warnMaxSizeExceeded ? "danger" : "none"}
              >
                <AnchorButton
                  disabled={!haveBothCellSets || warnMaxSizeExceeded}
                  intent={warnMaxSizeExceeded ? "danger" : "primary"}
                  data-testid="diffexp-button"
                  loading={differential.loading}
                  icon="left-join"
                  fill
                  onClick={this.computeDiffExp}
                />
              </Tooltip>
            }
            content={
              <div
                style={{
                  display: "flex",
                  justifyContent: "flex-start",
                  alignItems: "flex-end",
                  flexDirection: "column",
                  padding: 10,
                  maxWidth: 310,
                }}
              >
                <p>
                  {`The total number of cells for differential expression computation
                may not exceed ${diffexpCellcountMax}`}
                </p>
                <Button
                  type="button"
                  data-testid="diffexp-maxsize-exceeded-warning-dismiss"
                  intent="warning"
                  onClick={this.clearDifferentialExpression}
                >
                  Dismiss and clear cell sets
                </Button>
                <Button
                  type="button"
                  data-testid="diffexp-popover-dismiss"
                  intent="none"
                  onClick={this.handlePopoverDismiss}
                >
                  Dismiss
                </Button>
              </div>
            }
          />
        ) : null}

        {differential.diffExp ? (
          <Tooltip
            content="Remove differentially expressed gene list and clear cell selections"
            position="bottom"
            hoverOpenDelay={globals.tooltipHoverOpenDelayQuick}
          >
            <Button
              type="button"
              fill
              intent="warning"
              onClick={this.clearDifferentialExpression}
            >
              Clear Differential Expression
            </Button>
          </Tooltip>
        ) : null}
      </ButtonGroup>
    );
  }
}

export default DiffexpButtons;
