import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import {
  Button,
  Tooltip,
  InputGroup,
  Menu,
  Dialog,
  MenuItem,
  Popover,
  Classes,
  Icon,
  Position,
  PopoverInteractionKind,
  Colors
} from "@blueprintjs/core";
import { Suggest, Select } from "@blueprintjs/select";
import fuzzysort from "fuzzysort";
import actions from "../../actions";

import * as globals from "../../globals";
import Value from "./value";
import sortedCategoryValues from "./util";
import { AnnotationsHelpers } from "../../util/stateManager";

const filterOntology = (query, genes) =>
  /* fires on load, once, and then for each character typed into the input */
  fuzzysort.go(query, genes, {
    limit: 5,
    threshold: -10000 // don't return bad results
  });

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
      isCreatingNewCategory,
      allCategoryNames
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
                intent={
                  this.categoryNameError(newCategoryText) ? "warning" : "none"
                }
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
                {this.categoryNameErrorMessage(newCategoryText)}
              </p>
            </div>

            <p>
              Optionally duplicate all labels & cell assignments from existing
              category into new category:
            </p>
            {this.children}
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
