import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../annoDialog";
import LabelInput from "../../labelInput";
import { labelPrompt } from "../labelUtil";

import { AnnotationsHelpers } from "../../../util/stateManager";
import actions from "../../../actions";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  annotations: (state as any).annotations,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  schema: (state as any).annoMatrix?.schema,
}))
// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
class AnnoDialogEditCategoryName extends React.PureComponent<{}, State> {
  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);
    this.state = {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type '{... Remove this comment to see the full error message
      newCategoryText: props.metadataField,
    };
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleChangeOrSelect = (name: any) => {
    this.setState({
      newCategoryText: name,
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  disableEditCategoryMode = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "annotation: disable category edit mode",
    });
    this.setState({ newCategoryText: metadataField });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleEditCategory = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField } = this.props;
    const { newCategoryText } = this.state;
    /*
        test for uniqueness against *all* annotation names, not just the subset
        we render as categorical.
        */
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'schema' does not exist on type 'Readonly... Remove this comment to see the full error message
    const { schema } = this.props;
    const allCategoryNames = schema.annotations.obs.columns.map(
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (c: any) => c.name
    );
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  editedCategoryNameError = (name: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
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
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'schema' does not exist on type 'Readonly... Remove this comment to see the full error message
    const { schema } = this.props;
    const allCategoryNames = schema.annotations.obs.columns.map(
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (c: any) => c.name
    );
    const categoryNameAlreadyExists = allCategoryNames.indexOf(name) > -1;
    const sameName = name === metadataField;
    if (categoryNameAlreadyExists && !sameName) {
      return "duplicate";
    }
    /* otherwise, no error */
    return false;
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  instruction = (name: any) => {
    return labelPrompt(
      this.editedCategoryNameError(name),
      "New, unique category name",
      ":"
    );
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  allCategoryNames() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'schema' does not exist on type 'Readonly... Remove this comment to see the full error message
    const { schema } = this.props;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    return schema.annotations.obs.columns.map((c: any) => c.name);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const { newCategoryText } = this.state;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
    const { metadataField, annotations } = this.props;
    return (
      <>
        <AnnoDialog
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ isActive: any; inputProps: { "data-testid"... Remove this comment to see the full error message
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
              // @ts-expect-error ts-migrate(2322) FIXME: Type '{ label: any; labelSuggestions: null; onChan... Remove this comment to see the full error message
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
