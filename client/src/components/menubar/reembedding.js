import React from "react";
import { connect } from "react-redux";
import { AnchorButton, ButtonGroup, Tooltip } from "@blueprintjs/core";
import * as globals from "../../globals";
import actions from "../../actions";
import styles from "./menubar.css";

@connect((state) => ({
  reembedController: state.reembedController,
  annoMatrix: state.annoMatrix,
}))
class Reembedding extends React.PureComponent {
  render() {
    const { dispatch, annoMatrix, reembedController } = this.props;
    const loading = !!reembedController?.pendingFetch;
    const disabled = annoMatrix.nObs === annoMatrix.schema.dataframe.nObs;
    const tipContent = disabled
      ? "Subset cells first, then click to recompute UMAP embedding."
      : "Click to recompute UMAP embedding on the current cell subset.";

    return (
      <ButtonGroup className={styles.menubarButton}>
        <Tooltip
          content={tipContent}
          position="bottom"
          hoverOpenDelay={globals.tooltipHoverOpenDelay}
        >
          <AnchorButton
            icon="new-object"
            disabled={disabled}
            onClick={() => dispatch(actions.requestReembed())}
            loading={loading}
          />
        </Tooltip>
      </ButtonGroup>
    );
  }
}

export default Reembedding;
