import React from "react";
import { connect } from "react-redux";
import actions from "../../actions";
import FilenameDialog from "./filenameDialog";

@connect((state) => ({
  annotations: state.annotations,
  saveInProgress: state.autosave?.saveInProgress ?? false,
  error: state.autosave?.error,
  writableCategoriesEnabled: state.config?.parameters?.annotations ?? false,
  annoMatrix: state.annoMatrix,
  lastSavedAnnoMatrix: state.autosave?.lastSavedAnnoMatrix,
}))
class Autosave extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      timer: null,
    };
  }

  componentDidMount() {
    const { writableCategoriesEnabled } = this.props;

    let { timer } = this.state;
    if (timer) clearInterval(timer);
    if (writableCategoriesEnabled) {
      timer = setInterval(this.tick, 2500);
    } else {
      timer = null;
    }
    this.setState({ timer });
  }

  componentWillUnmount() {
    const { timer } = this.state;
    if (timer) this.clearInterval(timer);
  }

  tick = () => {
    const { dispatch } = this.props;
    if (this.needToSave()) {
      dispatch(actions.saveObsAnnotationsAction());
    }
  };

  needToSave = () => {
    /* return true if we need to save, false if we don't */
    const { annoMatrix, lastSavedAnnoMatrix } = this.props;
    return actions.needToSaveObsAnnotations(annoMatrix, lastSavedAnnoMatrix);
  };

  statusMessage() {
    const { error, saveInProgress } = this.props;

    if (saveInProgress) return `Saving...`;

    if (error) {
      return `Autosave error: ${error}`;
    }

    return this.needToSave() ? "Unsaved" : "All saved";
  }

  render() {
    const {
      writableCategoriesEnabled,
      saveInProgress,
      lastSavedAnnoMatrix,
    } = this.props;
    const initialDataLoadComplete = lastSavedAnnoMatrix;

    if (!writableCategoriesEnabled) return null;

    return (
      <div
        id="autosave"
        data-testclass={
          !initialDataLoadComplete
            ? "autosave-init"
            : this.needToSave() || saveInProgress
            ? "autosave-incomplete"
            : "autosave-complete"
        }
        style={{
          position: "absolute",
          display: "inherit",
          right: 8,
          bottom: 8,
          zIndex: 1,
        }}
      >
        {this.statusMessage()}
        <FilenameDialog />
      </div>
    );
  }
}

export default Autosave;
