import React from "react";
import { connect } from "react-redux";
import * as globals from "../../globals";
import actions from "../../actions";
import FilenameDialog from "./filenameDialog";

@connect((state) => ({
  universe: state.universe,
  annotations: state.annotations,
  obsAnnotations: state.universe.obsAnnotations,
  saveInProgress: state.autosave?.saveInProgress ?? false,
  lastSavedObsAnnotations: state.autosave?.lastSavedObsAnnotations,
  error: state.autosave?.error,
  writableCategoriesEnabled: state.config?.parameters?.["annotations"] ?? false,
  initialDataLoadComplete: state.autosave?.initialDataLoadComplete,
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
    const { dispatch, saveInProgress } = this.props;
    if (this.needToSave() && !saveInProgress) {
      dispatch(actions.saveObsAnnotations());
    }
  };

  needToSave = () => {
    /* return true if we need to save, false if we don't */
    const { obsAnnotations, lastSavedObsAnnotations } = this.props;
    return (
      lastSavedObsAnnotations && obsAnnotations !== lastSavedObsAnnotations
    );
  };

  statusMessage() {
    const { error } = this.props;
    if (error) {
      return `Autosave error: ${error}`;
    }
    return this.needToSave() ? "Unsaved" : "All saved";
  }

  render() {
    const {
      writableCategoriesEnabled,
      saveInProgress,
      initialDataLoadComplete,
    } = this.props;
    return writableCategoriesEnabled ? (
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
    ) : null;
  }
}

export default Autosave;
