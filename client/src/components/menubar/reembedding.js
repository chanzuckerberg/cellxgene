import React from "react";
import { connect } from "react-redux";
import {
  AnchorButton,
  ButtonGroup,
  Tooltip,
  Dialog,
  ControlGroup,
  Button
} from "@blueprintjs/core";
import * as globals from "../../globals";
import actions from "../../actions";
import styles from "./menubar.css";
import DimredPanel from "./dimredpanel";

@connect((state) => ({
  reembedController: state.reembedController,
  preprocessController: state.preprocessController,
  reembedParams: state.reembedParameters,
  annoMatrix: state.annoMatrix,
  idhash: state.config?.parameters?.["annotations-user-data-idhash"] ?? null,
  obsCrossfilter: state.obsCrossfilter,
  layoutChoice: state.layoutChoice
}))
class Reembedding extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      setReembedDialogActive: false,
      embName: "",
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
    const { dispatch, reembedParams, layoutChoice, obsCrossfilter } = this.props;
    const { embName } = this.state
    let parentName;
    if (obsCrossfilter.countSelected() === obsCrossfilter.annoMatrix.nObs) {
      if (layoutChoice.current.includes(";;")){
        parentName = layoutChoice.current.split(";;")
        parentName.pop()
        parentName = parentName.join(';;');
      } else{
        parentName="";
      }
    } else {
      parentName = layoutChoice.current;
    }
    dispatch(actions.requestReembed(reembedParams,parentName, embName));
    this.setState({
      setReembedDialogActive: false,
      embName: ""
    });
    // this is where you need to trigger subset if cells were filtered.
  };
  onNameChange = (name) => {
    this.setState({embName: name.target.value})
  }
  render() {
    const { setReembedDialogActive, embName } = this.state;
    const { reembedController, idhash, annoMatrix, obsCrossfilter, preprocessController, reembedParams } = this.props;
    const loading = !!reembedController?.pendingFetch || !!preprocessController?.pendingFetch;
    const tipContent =
      "Click to recompute UMAP embedding on the currently selected cells.";

    return (
      <ButtonGroup className={styles.menubarButton}>
        <Dialog
          icon="info-sign"
          onClose={this.handleDisableReembedDialog}
          title={`Reembedding on ${obsCrossfilter.countSelected()}/${annoMatrix.schema.dataframe.nObs} cells.`}
          autoFocus
          canEscapeKeyClose
          canOutsideClickClose
          enforceFocus
          initialStepIndex={0}
          isOpen={setReembedDialogActive}
          usePortal
        >        
          <div style={{
            marginLeft: "10px",
            marginRight: "10px"
          }}>        
            <DimredPanel embName={embName} onChange={this.onNameChange} idhash={idhash} />
            <ControlGroup style={{paddingTop: "15px"}} fill={true} vertical={false}>
                <Button onClick={this.handleDisableReembedDialog}>Close</Button>
                <Button disabled={reembedParams.doBatch && reembedParams.batchKey===""} onClick={this.handleRunAndDisableReembedDialog} intent="primary"> Run </Button>                 
            </ControlGroup>            
          </div>
        </Dialog>
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
