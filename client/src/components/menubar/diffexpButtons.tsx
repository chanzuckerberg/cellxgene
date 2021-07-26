import React from "react";
import { connect } from "react-redux";
import { ButtonGroup, AnchorButton, Tooltip } from "@blueprintjs/core";
import * as globals from "../../globals";
// @ts-expect-error ts-migrate(2307) FIXME: Cannot find module './menubar.css' or its correspo... Remove this comment to see the full error message
import styles from "./menubar.css";
import actions from "../../actions";
import CellSetButton from "./cellSetButtons";

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  differential: (state as any).differential,
  diffexpMayBeSlow:
    (state as any).config?.parameters?.["diffexp-may-be-slow"] ?? false,
  diffexpCellcountMax: (state as any).config?.limits?.diffexp_cellcount_max,
}))
class DiffexpButtons extends React.PureComponent {
  computeDiffExp = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
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
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'differential' does not exist on type 'Re... Remove this comment to see the full error message
    const { differential, diffexpMayBeSlow, diffexpCellcountMax } = this.props;
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
      <ButtonGroup className={styles.menubarButton}>
        {/* @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call. */}
        <CellSetButton eitherCellSetOneOrTwo={1} />
        {/* @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call. */}
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
