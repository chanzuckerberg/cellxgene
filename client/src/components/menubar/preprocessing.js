import React from "react";
import { connect } from "react-redux";
import {
  AnchorButton,
  ButtonGroup,
  Button,
  ControlGroup,
  Tooltip,
  Dialog,
} from "@blueprintjs/core";
import * as globals from "../../globals";
import actions from "../../actions";
import styles from "./menubar.css";
import PrepPanel from "./preppanel";

@connect((state) => ({
  reembedController: state.reembedController,
  preprocessController: state.preprocessController,
  reembedParams: state.reembedParameters,
  idhash: state.config?.parameters?.["annotations-user-data-idhash"] ?? null,
}))
class Preprocessing extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      setPreprocessingDialogActive: false,
    };
  }

  handleEnablePreprocessingDialog = () => {
    this.setState({ setPreprocessingDialogActive: true });
  };

  handleDisablePreprocessingDialog = () => {
    this.setState({
      setPreprocessingDialogActive: false,
    });
  };

  handleRunAndDisablePreprocessingDialog = () => {
    const { dispatch, reembedParams } = this.props;
    
    dispatch(actions.requestPreprocessing(reembedParams));
    this.setState({
      setPreprocessingDialogActive: false,
    });
    // this is where you need to trigger subset if cells were filtered.
  };
  onNameChange = (name) => {
    this.setState({embName: name.target.value})
  }
  render() {
    const { setPreprocessingDialogActive } = this.state;
    const { reembedController, preprocessController, idhash } = this.props;
    const loading = !!reembedController?.pendingFetch || !!preprocessController?.pendingFetch;    
    const tipContent =
      "Click to perform preprocessing.";
    return (
      <ButtonGroup className={styles.menubarButton}>
        <Dialog
          icon="info-sign"
          onClose={this.handleDisablePreprocessingDialog}
          title={`Preprocessing`}
          autoFocus
          canEscapeKeyClose
          canOutsideClickClose
          enforceFocus
          initialStepIndex={0}
          isOpen={setPreprocessingDialogActive}
          usePortal
        >
          <div style={{
            marginLeft: "10px",
            marginRight: "10px"
          }}>
            <PrepPanel idhash={idhash} />
            <ControlGroup style={{paddingTop: "15px"}} fill={true} vertical={false}>
              <Button onClick={this.handleDisablePreprocessingDialog}>Close</Button>
              <Button onClick={this.handleRunAndDisablePreprocessingDialog} intent="primary"> Preprocess </Button>                 
            </ControlGroup>            
          </div>
        </Dialog>
        <Tooltip
          content={tipContent}
          position="bottom"
          hoverOpenDelay={globals.tooltipHoverOpenDelay}
        >
          <AnchorButton
            icon="clean"
            loading={loading}
            onClick={this.handleEnablePreprocessingDialog}
            data-testid="preprocessing-options"
          />
        </Tooltip>
      </ButtonGroup>
    );
  }
}

export default Preprocessing;
