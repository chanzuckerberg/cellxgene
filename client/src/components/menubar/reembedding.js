import React from "react";
import { connect } from "react-redux";
import {
  AnchorButton,
  ButtonGroup,
  Tooltip,
  DialogStep,
  MultistepDialog,
} from "@blueprintjs/core";
import * as globals from "../../globals";
import actions from "../../actions";
import styles from "./menubar.css";
import PrepPanel from "./preppanel";
import BatchPanel from "./batchpanel";
import DimredPanel from "./dimredpanel";

@connect((state) => ({
  reembedController: state.reembedController,
  reembedParams: state.reembedParameters,
  annoMatrix: state.annoMatrix,
  idhash: state.config?.parameters?.["annotations-user-data-idhash"] ?? null,
}))
class Reembedding extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      setReembedDialogActive: false,
    };
  }

  handleEnableReembedDialog = () => {
    this.setState({ setReembedDialogActive: true });
  };

  handleDisableReembedDialog = () => {
    this.setState({
      setReembedDialogActive: false,
    });
  };

  handleRunAndDisableReembedDialog = () => {
    const { dispatch, reembedParams } = this.props;
    dispatch(actions.requestReembed(reembedParams));
    this.setState({
      setReembedDialogActive: false,
    });
    // this is where you need to trigger subset if cells were filtered.
  };

  render() {
    const { setReembedDialogActive } = this.state;

    const { reembedController, idhash, reembedParams, annoMatrix } = this.props;
    const loading = !!reembedController?.pendingFetch;
    const tipContent =
      "Click to recompute UMAP embedding on the currently selected cells.";
    const finalButtonProps = {
      intent: "primary",
      onClick: this.handleRunAndDisableReembedDialog,
      text: "Run",
    };
    const nextButtonProps = {
      intent: "primary",
      disabled: reembedParams.doBatch && reembedParams.batchKey === "",
      text: "Next",
    };
    return (
      <ButtonGroup className={styles.menubarButton}>
        <MultistepDialog
          icon="info-sign"
          onClose={this.handleDisableReembedDialog}
          finalButtonProps={finalButtonProps}
          title={`Reembedding on ${annoMatrix.nObs}/${annoMatrix.schema.dataframe.nObs} cells.`}
          autoFocus
          canEscapeKeyClose
          canOutsideClickClose
          enforceFocus
          initialStepIndex={0}
          isOpen={setReembedDialogActive}
          usePortal
        >
          <DialogStep
            id="preprocessing"
            panel={<PrepPanel idhash={idhash} />}
            title="Preprocessing"
          />
          <DialogStep
            id="batchcorrect"
            panel={
              <div>
                <BatchPanel idhash={idhash} />
              </div>
            }
            title="Batch correction"
            nextButtonProps={nextButtonProps}
          />
          <DialogStep
            id="dimred"
            panel={
              <div>
                <DimredPanel idhash={idhash} />
              </div>
            }
            title="Dimensionality reduction"
          />
        </MultistepDialog>
        <Tooltip
          content={tipContent}
          position="bottom"
          hoverOpenDelay={globals.tooltipHoverOpenDelay}
        >
          <AnchorButton
            icon="new-object"
            loading={loading}
            onClick={this.handleEnableReembedDialog}
            data-testid="reembedding-options"
          />
        </Tooltip>
      </ButtonGroup>
    );
  }
}

export default Reembedding;
