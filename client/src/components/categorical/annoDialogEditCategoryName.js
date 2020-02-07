import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import { Colors } from "@blueprintjs/core";
import AnnoDialog from "./annoDialog";
import AnnoInputs from "./annoInputs";

import { AnnotationsHelpers } from "../../util/stateManager";

@connect(state => ({
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  universe: state.universe,
  ontology: state.ontology,
  ontologyLoading: state.ontology?.loading,
  ontologyEnabled: state.ontology?.enabled
}))
class AnnoDialogEditCategoryName extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      newCategoryText: props.metadataField
    };
  }

  handleCategoryEditTextChange = txt => {
    this.setState({
      newCategoryText: txt
    });
  };

  disableEditCategoryMode = () => {
    const { dispatch } = this.props;
    dispatch({
      type: "annotation: disable category edit mode"
    });
  };

  handleEditCategory = () => {
    const { dispatch, metadataField, categoricalSelection } = this.props;
    const { newCategoryText } = this.state;

    const allCategoryNames = _.keys(categoricalSelection);

    if (
      (allCategoryNames.indexOf(newCategoryText) > -1 &&
        newCategoryText !== metadataField) ||
      newCategoryText === ""
    ) {
      return;
    }

    dispatch({
      type: "annotation: category edited",
      metadataField,
      newCategoryText,
      data: newCategoryText
    });
  };

  categoryNameErrorMessage = () => {
    const err = this.editedCategoryNameError();
    if (err === false) return null;

    const errorMessageMap = {
      /* map error code to human readable error message */
      "empty-string": "Blank names not allowed",
      duplicate: "Category name must be unique",
      "trim-spaces": "Leading and trailing spaces not allowed",
      "illegal-characters":
        "Only alphanumeric and special characters (-_.) allowed",
      "multi-space-run": "Multiple consecutive spaces not allowed"
    };
    const errorMessage = errorMessageMap[err] ?? "error";
    return (
      <span
        style={{
          display: "block",
          fontStyle: "italic",
          fontSize: 12,
          marginTop: 5,
          color: Colors.ORANGE3
        }}
      >
        {errorMessage}
      </span>
    );
  };

  editedCategoryNameError = () => {
    const { metadataField, categoricalSelection } = this.props;
    const { newCategoryText } = this.state;

    /* check for syntax errors in category name */
    const error = AnnotationsHelpers.annotationNameIsErroneous(newCategoryText);
    if (error) {
      return error;
    }

    /* check for duplicative categories */
    const allCategoryNames = _.keys(categoricalSelection);
    const categoryNameAlreadyExists =
      allCategoryNames.indexOf(newCategoryText) > -1;
    const sameName = newCategoryText === metadataField;
    if (categoryNameAlreadyExists && !sameName) {
      return "duplicate";
    }

    /* otherwise, no error */
    return false;
  };

  render() {
    const { newCategoryText } = this.state;
    const { metadataField, annotations } = this.props;

    return (
      <>
        <AnnoDialog
          isActive={
            annotations.isEditingCategoryName &&
            annotations.categoryBeingEdited === metadataField
          }
          inputProps={{ "data-testid": `${metadataField}:edit-category-name-dialog` }}
          primaryButtonProps={{ "data-testid": `${metadataField}:submit-category-edit` }}
          title="Edit category name"
          instruction="New, unique category name:"
          cancelTooltipContent="Close this dialog without editing this category."
          primaryButtonText="Edit category name"
          text={newCategoryText}
          validationError={this.editedCategoryNameError(newCategoryText)}
          errorMessage={this.categoryNameErrorMessage(newCategoryText)}
          handleSubmit={this.handleEditCategory}
          handleCancel={this.disableEditCategoryMode}
          annoInput={
            <AnnoInputs
              inputProps={{ "data-testid": `${metadataField}:edit-category-name-text`}}
              useSuggest={false}
              text={newCategoryText}
              handleTextChange={this.handleCategoryEditTextChange}
            />
          }
        />
      </>
    );
  }
}

export default AnnoDialogEditCategoryName;
