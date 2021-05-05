import React from "react";
import { connect } from "react-redux";

import {
  Button,
  Classes,
  Code,
  Colors,
  Dialog,
  InputGroup,
  Tooltip,
} from "@blueprintjs/core";

@connect((state) => ({
  idhash: state.config?.parameters?.["annotations-user-data-idhash"] ?? null,
  annotations: state.annotations,
  auth: state.config?.authentication,
  userInfo: state.userInfo,
  writableCategoriesEnabled: state.config?.parameters?.annotations ?? false,
  writableGenesetsEnabled: !(
    state.config?.parameters?.annotations_genesets_readonly ?? true
  ),
}))
class FilenameDialog extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      filenameText: "",
    };
  }

  dismissFilenameDialog = () => {};

  handleCreateFilename = () => {
    const { dispatch } = this.props;
    const { filenameText } = this.state;

    dispatch({
      type: "set annotations collection name",
      data: filenameText,
    });
  };

  filenameError = () => {
    const legalNames = /^\w+$/;
    const { filenameText } = this.state;
    let err = false;

    if (filenameText === "") {
      err = "empty_string";
    } else if (!legalNames.test(filenameText)) {
      /*
      IMPORTANT: this test must ultimately match the test applied by the
      backend, which is designed to ensure a safe file name can be created
      from the data collection name.  If you change this, you will also need
      to change the validation code in the backend, or it will have no effect.
      */
      err = "characters";
    }

    return err;
  };

  filenameErrorMessage = () => {
    const err = this.filenameError();
    let markup = null;

    if (err === "empty_string") {
      markup = (
        <span
          style={{
            fontStyle: "italic",
            fontSize: 12,
            marginTop: 5,
            color: Colors.ORANGE3,
          }}
        >
          Name cannot be blank
        </span>
      );
    } else if (err === "characters") {
      markup = (
        <span
          style={{
            fontStyle: "italic",
            fontSize: 12,
            marginTop: 5,
            color: Colors.ORANGE3,
          }}
        >
          Only alphanumeric and underscore allowed
        </span>
      );
    }
    return markup;
  };

  render() {
    const {
      writableCategoriesEnabled,
      writableGenesetsEnabled,
      annotations,
      idhash,
      userInfo,
    } = this.props;
    const { filenameText } = this.state;

    return (writableCategoriesEnabled || writableGenesetsEnabled) &&
      annotations.promptForFilename &&
      !annotations.dataCollectionNameIsReadOnly &&
      !annotations.dataCollectionName &&
      userInfo.is_authenticated ? (
      <Dialog
        icon="tag"
        title="User Generated Data Directory"
        isOpen={!annotations.dataCollectionName}
        onClose={this.dismissFilenameDialog}
      >
        <form
          onSubmit={(e) => {
            e.preventDefault();
            this.handleCreateFilename();
          }}
        >
          <div className={Classes.DIALOG_BODY} data-testid="annotation-dialog">
            <div style={{ marginBottom: 20 }}>
              <p>Name your user generated data directory:</p>
              <InputGroup
                autoFocus
                value={filenameText}
                intent={this.filenameError(filenameText) ? "warning" : "none"}
                onChange={(e) =>
                  this.setState({ filenameText: e.target.value })
                }
                leftIcon="tag"
                data-testid="new-annotation-name"
              />
              <p
                style={{
                  marginTop: 7,
                  visibility: this.filenameError(filenameText)
                    ? "visible"
                    : "hidden",
                  color: Colors.ORANGE3,
                }}
              >
                {this.filenameErrorMessage(filenameText)}
              </p>
            </div>
            <div>
              <p>
                Your annotations are stored in this file:
                <Code>
                  {filenameText}-{idhash}.csv
                </Code>
              </p>
              <p>
                Your gene sets are stored in this file:
                <Code>
                  {filenameText}-gene-sets-{idhash}.csv
                </Code>
              </p>
              <p style={{ fontStyle: "italic" }}>
                (We added a unique ID to your filename)
              </p>
            </div>
          </div>
          <div className={Classes.DIALOG_FOOTER}>
            <div className={Classes.DIALOG_FOOTER_ACTIONS}>
              <Tooltip content="Cancel naming collection">
                <Button onClick={this.dismissFilenameDialog}>Cancel</Button>
              </Tooltip>
              <Button
                disabled={!filenameText || this.filenameError(filenameText)}
                onClick={this.handleCreateFilename}
                intent="primary"
                type="submit"
                data-testid="submit-annotation"
              >
                Create user generated data directory
              </Button>
            </div>
          </div>
        </form>
      </Dialog>
    ) : null;
  }
}

export default FilenameDialog;
