import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../annoDialog";
import LabelInput from "../../labelInput";
import { labelPrompt } from "../labelUtil";

import { AnnotationsHelpers } from "../../../util/stateManager";
import actions from "../../../actions";

@connect((state) => ({
  annotations: state.annotations,
  schema: state.annoMatrix?.schema,
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
    const { dispatch, metadataField } = this.props;
    const { newCategoryText } = this.state;

    /*
    test for uniqueness against *all* annotation names, not just the subset
    we render as categorical.
    */
    const { schema } = this.props;
    const allCategoryNames = schema.annotations.obs.columns.map((c) => c.name);

    if (
      (allCategoryNames.indexOf(newCategoryText) > -1 &&
        newCategoryText !== metadataField) ||
      newCategoryText === ""
    ) {
      return;
    }

    this.disableEditCategoryMode();

    if (metadataField !== newCategoryText)
      dispatch(
        actions.annotationRenameCategoryAction(metadataField, newCategoryText)
      );
    e.preventDefault();
  };

  editedCategoryNameError = (name) => {
    const { metadataField } = this.props;

    /* check for syntax errors in category name */
    const error = AnnotationsHelpers.annotationNameIsErroneous(name);
    if (error) {
      return error;
    }

    /* check for duplicative categories */

    /*
    test for uniqueness against *all* annotation names, not just the subset
    we render as categorical.
    */
    const { schema } = this.props;
    const allCategoryNames = schema.annotations.obs.columns.map((c) => c.name);

    const categoryNameAlreadyExists = allCategoryNames.indexOf(name) > -1;
    const sameName = name === metadataField;
    if (categoryNameAlreadyExists && !sameName) {
      return "duplicate";
    }

    /* otherwise, no error */
    return false;
  };

  instruction = (name) => labelPrompt(
      this.editedCategoryNameError(name),
      "New, unique category name",
      ":"
    );

  allCategoryNames() {
    const { schema } = this.props;
    return schema.annotations.obs.columns.map((c) => c.name);
  }

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
              labelSuggestions={null}
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
