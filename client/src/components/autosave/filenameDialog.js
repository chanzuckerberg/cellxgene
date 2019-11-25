import React from "react";
import { connect } from "react-redux";

import {
  Button,
  Tooltip,
  InputGroup,
  Dialog,
  Classes,
  Colors
} from "@blueprintjs/core";

@connect(state => ({
  universe: state.universe,
  annotations: state.annotations,
  obsAnnotations: state.universe.obsAnnotations,
  saveInProgress: state.autosave?.saveInProgress ?? false,
  lastSavedObsAnnotations: state.autosave?.lastSavedObsAnnotations,
  error: state.autosave?.error,
  writableCategoriesEnabled: state.config?.parameters?.["annotations"] ?? false
}))
class FilenameDialog extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      filenameText: ""
    };
  }

  dismissFilenameDialog = () => {};

  handleCreateFilename = () => {
    const { dispatch } = this.props;
    const { filenameText } = this.state;

    dispatch({
      type: "set annotations filename",
      data: filenameText
    });
  };

  filenameError = () => {
    const { filenameText } = this.state;
    let err;

    if (filenameText === "") {
      err = "empty_string";
    }

    /* add legal character check */

    return err;
  };

  filenameErrorMessage = () => {
    const err = this.filenameError();
    if (!err) return null;

    let markup = null;

    if (err === "empty_string") {
      markup = (
        <span
          style={{
            fontStyle: "italic",
            fontSize: 12,
            marginTop: 5,
            color: Colors.ORANGE3
          }}
        >
          {"Filename cannot be blank"}
        </span>
      );
    } else if (err === "characters") {
      markup = (
        <span
          style={{
            fontStyle: "italic",
            fontSize: 12,
            marginTop: 5,
            color: Colors.ORANGE3
          }}
        >
          {"Only alphanumeric and underscore allowed"}
        </span>
      );
    }
    return markup;
  };

  render() {
    const { writableCategoriesEnabled, annotations } = this.props;
    const { filenameText } = this.state;

    return writableCategoriesEnabled &&
      !annotations.dataCollectionNameIsReadOnly &&
      !annotations.dataCollectionName ? (
      <Dialog
        icon="tag"
        title="Filename"
        isOpen={!annotations.dataCollectionName}
        onClose={this.dismissFilenameDialog}
      >
        <form
          onSubmit={e => {
            e.preventDefault();
            this.handleCreateFilename();
          }}
        >
          <div className={Classes.DIALOG_BODY}>
            <div style={{ marginBottom: 20 }}>
              <p>Set a filename:</p>
              <InputGroup
                autoFocus
                value={filenameText}
                intent={this.filenameError(filenameText) ? "warning" : "none"}
                onChange={e => this.setState({ filenameText: e.target.value })}
                leftIcon="tag"
              />
              <p
                style={{
                  marginTop: 7,
                  visibility: this.filenameError(filenameText)
                    ? "visible"
                    : "hidden",
                  color: Colors.ORANGE3
                }}
              >
                {this.filenameErrorMessage(filenameText)}
              </p>
            </div>
          </div>
          <div className={Classes.DIALOG_FOOTER}>
            <div className={Classes.DIALOG_FOOTER_ACTIONS}>
              <Tooltip content="Cancel naming CSV">
                <Button onClick={this.dismissFilenameDialog}>Cancel</Button>
              </Tooltip>
              <Button
                disabled={!filenameText || this.filenameError(filenameText)}
                onClick={this.handleCreateFilename}
                intent="primary"
                type="submit"
              >
                Create CSV
              </Button>
            </div>
          </div>
        </form>
      </Dialog>
    ) : null;
  }
}

export default FilenameDialog;
