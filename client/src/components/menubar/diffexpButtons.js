// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import { Popover, Button, ButtonGroup, AnchorButton, Tooltip, Position } from "@blueprintjs/core";
import * as globals from "../../globals";
import actions from "../../actions";
import CellSetButton from "./cellSetButtons";

@connect(state => ({
  config: state.config,
  crossfilter: state.crossfilter,
  differential: state.differential,
  celllist1: state.differential?.celllist1,
  celllist2: state.differential?.celllist2,
  diffexpMayBeSlow: state.config?.parameters?.["diffexp-may-be-slow"] ?? false,
  diffexpCellcountMax: state.config?.limits?.diffexp_cellcount_max
}))
class DiffexpButtons extends React.Component {
  static isValidDigitKeyEvent(e) {
    /*
    Return true if this event is necessary to enter a percent number input.
    Return false if not.

    Returns true for events with keys: backspace, control, alt, meta, [0-9],
    or events that don't have a key.
    */
    if (e.key === null) return true;
    if (e.ctrlKey || e.altKey || e.metaKey) return true;

    // concept borrowed from blueprint's numericInputUtils:
    // keys that print a single character when pressed have a `key` name of
    // length 1. every other key has a longer `key` name (e.g. "Backspace",
    // "ArrowUp", "Shift"). since none of those keys can print a character
    // to the field--and since they may have important native behaviors
    // beyond printing a character--we don't want to disable their effects.
    const isSingleCharKey = e.key.length === 1;
    if (!isSingleCharKey) return true;

    const key = e.key.charCodeAt(0) - 48; /* "0" */
    return key >= 0 && key <= 9;
  }

  constructor(props) {
    super(props);
    this.state = {
      pendingClipPercentiles: null
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
      diffExp: differential.diffExp
    });
    dispatch({
      type: "clear scatterplot"
    });
  };

  render() {

    /* diffexp-related buttons may be disabled */
    const { differential, diffexpMayBeSlow, diffexpCellcountMax } = this.props;

    const haveBothCellSets =
      !!differential.celllist1 && !!differential.celllist2;
    
    const haveEitherCellSet = !!differential.celllist1 || !!differential.celllist2;

    const slowMsg = diffexpMayBeSlow
      ? " (CAUTION: large dataset - may take longer or fail)"
      : "";
    const tipMessage = `See top 10 differentially expressed genes${slowMsg}`; 

    const warnMaxSizeExceeded = haveEitherCellSet &&/* diffexpCellcountMax &&*/
      (
        (differential.celllist1?.length || 0) + 
        (differential.celllist2?.length || 0) > 
        diffexpCellcountMax
      )
    
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
          isOpen={warnMaxSizeExceeded}
          position={Position.BOTTOM}
          target={
            <Tooltip
              content={tipMessage}
              position="bottom"
              hoverOpenDelay={globals.tooltipHoverOpenDelayQuick}
            >
              <AnchorButton
                disabled={!haveBothCellSets}
                intent="primary"
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
                maxWidth: 310
              }}>
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
    )
  }
}

export default DiffexpButtons;
