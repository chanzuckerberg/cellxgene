import React from "react";
import { connect } from "react-redux";
import { ButtonGroup, AnchorButton, Tooltip } from "@blueprintjs/core";
import * as globals from "../../globals";
import styles from "./menubar.css";
import actions from "../../actions";
import CellSetButton from "./cellSetButtons";

@connect((state) => ({
  differential: state.differential,
  diffexpMayBeSlow: state.config?.parameters?.["diffexp-may-be-slow"] ?? false,
  diffexpCellcountMax: state.config?.limits?.diffexp_cellcount_max,
}))
class DiffexpButtons extends React.PureComponent {
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

  render() {
    /* diffexp-related buttons may be disabled */
    const { differential, diffexpMayBeSlow, diffexpCellcountMax } = this.props;

    const haveBothCellSets =
      !!differential.celllist1 && !!differential.celllist2;

    const haveEitherCellSet =
      !!differential.celllist1 || !!differential.celllist2;

    const slowMsg = diffexpMayBeSlow
      ? " (CAUTION: large dataset - may take longer or fail)"
      : "";
    const tipMessage = `See top differentially expressed genes${slowMsg}`;
    const tipMessageWarn = `The total number of cells for differential expression computation
                            may not exceed ${diffexpCellcountMax}. Try reselecting new protein sets.`;

    const warnMaxSizeExceeded =
      haveEitherCellSet &&
      !!diffexpCellcountMax &&
      (differential.celllist1?.length ?? 0) +
        (differential.celllist2?.length ?? 0) >
        diffexpCellcountMax;

    return (
      <ButtonGroup className={styles.menubarButton}>
        <CellSetButton eitherCellSetOneOrTwo={1} />
        <CellSetButton eitherCellSetOneOrTwo={2} />
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
      </ButtonGroup>
    );
  }
}

export default DiffexpButtons;
