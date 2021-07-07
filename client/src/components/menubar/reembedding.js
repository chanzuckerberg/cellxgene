import React from "react";
import { connect } from "react-redux";
import {
  AnchorButton,
  ButtonGroup,
  Tooltip,
  Slider,
  MenuItem,
  DialogStep,
  MultistepDialog,
  Collapse
} from "@blueprintjs/core";
import * as globals from "../../globals";
import actions from "../../actions";
import styles from "./menubar.css";
import PrepPanel from "./preppanel";
import BatchPanel from "./batchpanel";
import DimredPanel from "./dimredpanel";

@connect((state) => ({
  reembedController: state.reembedController,
  annoMatrix: state.annoMatrix,
  dimredParams: state.reembedParameters.dimredParams,
  prepParams: state.reembedParameters.prepParams,
  batchCorrectionParams: state.reembedParameters.batchParams,
}))
class Reembedding extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      setReembedDialogActive: false,
      params: {}
    };
  }  
  handleEnableReembedDialog = () => {
    this.setState({ setReembedDialogActive: true });
  };


  handleDisableReembedDialog = () => {
    this.setState({
      setReembedDialogActive: false
    });
  };

  handleRunAndDisableReembedDialog = () => {
    const { dispatch } = this.props;
    dispatch(actions.requestReembed(this.state.params));
    this.setState({
      setReembedDialogActive: false
    });
  };
  render() {
    const {
      setReembedDialogActive, params
    } = this.state;

    const { dispatch, annoMatrix, reembedController, dimredParams, prepParams, batchCorrectionParams } = this.props;
    const loading = !!reembedController?.pendingFetch;
    const tipContent = "Click to recompute UMAP embedding on the currently selected cells.";
      const finalButtonProps = {
        intent: "primary",
        onClick: this.handleRunAndDisableReembedDialog,
        text: "Run",
    };
    const dialogparams = {
      autoFocus: true,
      canEscapeKeyClose: true,
      canOutsideClickClose: true,
      enforceFocus: true,
      initialStepIndex: 0,
      isOpen: setReembedDialogActive,
      usePortal: true,
    };
    return (
      <ButtonGroup className={styles.menubarButton}>
        <MultistepDialog
            icon="info-sign"
            onClose={this.handleDisableReembedDialog}
            finalButtonProps={finalButtonProps}
            title="Reembedding parameters"
            {...dialogparams}
        >
            <DialogStep
                id="preprocessing"
                panel={
                  <PrepPanel/>
                }
                title="Preprocessing"
            />            
            <DialogStep
                id="batchcorrect"
                panel={
                  <div>
                    <BatchPanel/>
                  </div>
                }
                title="Batch correction"
            />
            <DialogStep
                id="dimred"
                panel={
                  <div>
                    <DimredPanel/>
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
