import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import AnnoDialog from "../annoDialog";
import LabelInput from "../labelInput";
import { labelPrompt } from "../labelUtil";

import { AnnotationsHelpers } from "../../../util/stateManager";

@connect((state) => ({
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  universe: state.universe,
  ontology: state.ontology,
}))
class AnnoDialogEditCategoryName extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      newCategoryText: props.metadataField,
    };
  }

  handleChangeOrSelect = (name) => {
    this.setState({
      newCategoryText: name,
    });
  };

  disableEditCategoryMode = () => {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "annotation: disable category edit mode",
    });
    this.setState({ newCategoryText: metadataField });
  };

  handleEditCategory = (e) => {
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

    this.disableEditCategoryMode();
    dispatch({
      type: "annotation: category edited",
      metadataField,
      newCategoryText,
      data: newCategoryText,
    });
    e.preventDefault();
  };

  editedCategoryNameError = (name) => {
    const { metadataField, categoricalSelection } = this.props;

    /* check for syntax errors in category name */
    const error = AnnotationsHelpers.annotationNameIsErroneous(name);
    if (error) {
      return error;
    }

    /* check for duplicative categories */
    const allCategoryNames = _.keys(categoricalSelection);
    const categoryNameAlreadyExists = allCategoryNames.indexOf(name) > -1;
    const sameName = name === metadataField;
    if (categoryNameAlreadyExists && !sameName) {
      return "duplicate";
    }

    /* otherwise, no error */
    return false;
  };

  instruction = (name) => {
    return labelPrompt(
      this.editedCategoryNameError(name),
      "New, unique category name",
      ":"
    );
  };

  render() {
    const { newCategoryText } = this.state;
    const { metadataField, annotations, ontology } = this.props;
    const ontologyEnabled = ontology?.enabled ?? false;

    return (
      <>
        <AnnoDialog
          isActive={
            annotations.isEditingCategoryName &&
            annotations.categoryBeingEdited === metadataField
          }
          inputProps={{
            "data-testid": `${metadataField}:edit-category-name-dialog`,
          }}
          primaryButtonProps={{
            "data-testid": `${metadataField}:submit-category-edit`,
          }}
          title="Edit category name"
          instruction={this.instruction(newCategoryText)}
          cancelTooltipContent="Close this dialog without editing this category."
          primaryButtonText="Edit category name"
          text={newCategoryText}
          validationError={this.editedCategoryNameError(newCategoryText)}
          handleSubmit={this.handleEditCategory}
          handleCancel={this.disableEditCategoryMode}
          annoInput={
            <LabelInput
              label={newCategoryText}
              labelSuggestions={ontologyEnabled ? ontology.terms : null}
              onChange={this.handleChangeOrSelect}
              onSelect={this.handleChangeOrSelect}
              inputProps={{
                "data-testid": `${metadataField}:edit-category-name-text`,
                leftIcon: "tag",
                intent: "none",
                autoFocus: true,
              }}
              newLabelMessage="New category"
            />
          }
        />
      </>
    );
  }
}

export default AnnoDialogEditCategoryName;
