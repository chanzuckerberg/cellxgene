// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { Button } from "@blueprintjs/core";
import { connect } from "react-redux";
import * as globals from "../../globals";
import Category from "./category";
import { AnnotationsHelpers } from "../../util/stateManager";
import AnnoDialog from "./annoDialog";
import AnnoInputs from "./annoInputs";
import AnnoSelect from "./annoSelect";

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
      categoryToDuplicate: null,
      activeSuggestItem: null
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

    /* 
    test for uniqueness against *all* annotation names, not just the subset
    we render as categorical.
    */
    const { schema } = this.props;
    const allCategoryNames = schema.annotations.obs.columns.map(c => c.name);

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
      "illegal-characters": "Only alphanumeric and special characters (-_.) allowed",
      "multi-space-run": "Multiple consecutive spaces not allowed"
    };
    const errorMessage = errorMessageMap[err] ?? "error";
    return <span>{errorMessage}</span>;
  };

  handleNewCategoryText = e => {
    this.setState({ newCategoryText: e.target.value });
  };

  handleChoice = e => {
    /* Blueprint Suggest format */
    this.setState({ newCategoryText: e.target });
  };

  handleSuggestActiveItemChange = item => {
    this.setState({ activeSuggestItem: item });
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
        <AnnoDialog
          isActive={createAnnoModeActive}
          title="Create new category"
          instruction="New, unique category name:"
          cancelTooltipContent="Close this dialog without creating a category."
          primaryButtonText="Create new category"
          text={newCategoryText}
          validationError={this.categoryNameError(newCategoryText)}
          errorMessage={this.categoryNameErrorMessage(newCategoryText)}
          handleSubmit={this.handleCreateUserAnno}
          handleCancel={this.handleDisableAnnoMode}
          annoInput={
            <AnnoInputs
              text={newCategoryText}
              handleItemChange={this.handleSuggestActiveItemChange}
              handleChoice={this.handleChoice}
              handleTextChange={this.handleNewCategoryText}
            />
          }
          annoSelect={
            <AnnoSelect
              handleModalDuplicateCategorySelection={
                this.handleModalDuplicateCategorySelection
              }
              categoryToDuplicate={categoryToDuplicate}
              allCategoryNames={allCategoryNames}
            />
          }
        />

        {/* READ ONLY CATEGORICAL FIELDS */}
        {/* this is duplicative but flat, could be abstracted */}
        {_.map(
          allCategoryNames,
          catName =>
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
        {_.map(
          allCategoryNames,
          catName =>
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
