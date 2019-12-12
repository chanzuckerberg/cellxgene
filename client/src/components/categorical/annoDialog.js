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
  colorAccessor: state.colors.colorAccessor,
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  universe: state.universe,
  ontology: state.ontology,
  ontologyLoading: state.ontology?.loading
}))
class AnnoDialog extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const {
      isActive,
      handleUserTyping,
      newCategoryText,
      errorMessage,
      validationError,
      annoSelect
    } = this.props;
    return (
      <Dialog
        icon="tag"
        title="Create new category"
        isOpen={isActive}
        onClose={this.handleDisableAnnoMode}
      >
        <form
          onSubmit={e => {
            e.preventDefault();
          }}
        >
          <div className={Classes.DIALOG_BODY}>
            <div style={{ marginBottom: 20 }}>
              <p>New, unique category name:</p>
              <InputGroup
                autoFocus
                value={newCategoryText}
                intent={validationError ? "warning" : "none"}
                onChange={e => handleUserTyping(e)}
                leftIcon="tag"
              />

              <p
                style={{
                  marginTop: 7,
                  visibility: this.categoryNameError(newCategoryText)
                    ? "visible"
                    : "hidden",
                  color: Colors.ORANGE3
                }}
              >
                {errorMessage}
              </p>
            </div>
            {annoSelect || null}
          </div>
          <div className={Classes.DIALOG_FOOTER}>
            <div className={Classes.DIALOG_FOOTER_ACTIONS}>
              <Tooltip content="Close this dialog without creating a category.">
                <Button onClick={this.handleDisableAnnoMode}>Cancel</Button>
              </Tooltip>
              <Button
                onClick={this.handleCreateUserAnno}
                disabled={
                  !newCategoryText || this.categoryNameError(newCategoryText)
                }
                intent="primary"
                type="submit"
              >
                Create new category
              </Button>
            </div>
          </div>
        </form>
      </Dialog>
    );
  }
}

export default AnnoDialog;
