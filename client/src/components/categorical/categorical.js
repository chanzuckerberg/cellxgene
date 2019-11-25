// jshint esversion: 6
import React from "react";
import _ from "lodash";
import {
  Button,
  Tooltip,
  InputGroup,
  Dialog,
  Classes,
  MenuItem,
  Colors
} from "@blueprintjs/core";
import { Select } from "@blueprintjs/select";
import { connect } from "react-redux";
import * as globals from "../../globals";
import Category from "./category";
import { AnnotationsHelpers } from "../../util/stateManager";

@connect(state => ({
  categoricalSelection: state.categoricalSelection,
  writableCategoriesEnabled: state.config?.parameters?.["annotations"] ?? false,
  schema: state.world?.schema
}))
class Categories extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      createAnnoModeActive: false,
      newCategoryText: "",
      categoryToDuplicate: null
    };
  }

  handleCreateUserAnno = () => {
    const { dispatch } = this.props;
    const { newCategoryText, categoryToDuplicate } = this.state;
    dispatch({
      type: "annotation: create category",
      data: newCategoryText,
      categoryToDuplicate
    });
    this.setState({
      createAnnoModeActive: false,
      categoryToDuplicate: null,
      newCategoryText: ""
    });
  };

  handleEnableAnnoMode = () => {
    this.setState({ createAnnoModeActive: true });
  };

  handleDisableAnnoMode = () => {
    this.setState({
      createAnnoModeActive: false,
      categoryToDuplicate: null,
      newCategoryText: ""
    });
  };

  handleModalDuplicateCategorySelection = d => {
    this.setState({ categoryToDuplicate: d });
  };

  categoryNameError = name => {
    /*
    return false if this is a LEGAL/acceptable category name or NULL/empty string,
    or return an error type.
    */

    /* allow empty string */
    if (name === "") return false;

    const { categoricalSelection } = this.props;
    const allCategoryNames = Object.keys(categoricalSelection);

    /* check category name syntax */
    const error = AnnotationsHelpers.annotationNameIsErroneous(name);
    if (error) {
      return error;
    }

    /* disallow duplicates */
    if (allCategoryNames.indexOf(name) !== -1) {
      return "duplicate";
    }

    /* otherwise, no error */
    return false;
  };

  categoryNameErrorMessage = name => {
    const err = this.categoryNameError(name);
    if (err === false) return null;

    const errorMessageMap = {
      /* map error code to human readable error message */
      "empty-string": "Blank names not allowed",
      duplicate: "Name must be unique",
      "trim-spaces": "Leading and trailing spaces not allowed",
      "illegal-characters": "Only alphanumeric, underscore and period allowed",
      "multi-space-run": "Multiple consecutive spaces not allowed"
    };
    const errorMessage = errorMessageMap[err] ?? "error";
    return <span>{errorMessage}</span>;
  };

  render() {
    const {
      createAnnoModeActive,
      categoryToDuplicate,
      newCategoryText
    } = this.state;
    const {
      categoricalSelection,
      writableCategoriesEnabled,
      schema
    } = this.props;
    if (!categoricalSelection) return null;

    /* all names, sorted in display order.  Will be rendered in this order */
    const allCategoryNames = Object.keys(categoricalSelection).sort();

    return (
      <div
        style={{
          padding: globals.leftSidebarSectionPadding
        }}
      >
        {/* READ ONLY CATEGORICAL FIELDS */}
        {/* this is duplicative but flat, could be abstracted */}
        {_.map(allCategoryNames, catName =>
          !schema.annotations.obsByName[catName].writable ? (
            <Category
              key={catName}
              metadataField={catName}
              createAnnoModeActive={createAnnoModeActive}
              isUserAnno={false}
            />
          ) : null
        )}
        {/* WRITEABLE FIELDS */}
        {_.map(allCategoryNames, catName =>
          schema.annotations.obsByName[catName].writable ? (
            <Category
              key={catName}
              metadataField={catName}
              createAnnoModeActive={createAnnoModeActive}
              isUserAnno
            />
          ) : null
        )}
        {writableCategoriesEnabled ? (
          <div>
            <Dialog
              icon="tag"
              title="Create new category"
              isOpen={createAnnoModeActive}
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
                        this.categoryNameError(newCategoryText)
                          ? "warning"
                          : "none"
                      }
                      onChange={e =>
                        this.setState({ newCategoryText: e.target.value })
                      }
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
                    Optionally duplicate all labels & cell assignments from
                    existing category into new category:
                  </p>
                  <Select
                    items={allCategoryNames}
                    filterable={false}
                    itemRenderer={(d, { handleClick }) => {
                      return (
                        <MenuItem onClick={handleClick} key={d} text={d} />
                      );
                    }}
                    noResults={<MenuItem disabled text="No results." />}
                    onItemSelect={d => {
                      this.handleModalDuplicateCategorySelection(d);
                    }}
                  >
                    {/* children become the popover target; render value here */}
                    <Button
                      text={
                        categoryToDuplicate || "None (all cells 'unassigned')"
                      }
                      rightIcon="double-caret-vertical"
                    />
                  </Select>
                </div>
                <div className={Classes.DIALOG_FOOTER}>
                  <div className={Classes.DIALOG_FOOTER_ACTIONS}>
                    <Tooltip content="Close this dialog without creating a category.">
                      <Button onClick={this.handleDisableAnnoMode}>
                        Cancel
                      </Button>
                    </Tooltip>
                    <Button
                      onClick={this.handleCreateUserAnno}
                      disabled={
                        !newCategoryText ||
                        this.categoryNameError(newCategoryText)
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
            <Button onClick={this.handleEnableAnnoMode} intent="primary">
              Create new category
            </Button>
          </div>
        ) : null}
      </div>
    );
  }
}

export default Categories;
