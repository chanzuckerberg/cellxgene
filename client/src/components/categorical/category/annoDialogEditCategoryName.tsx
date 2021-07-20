import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../annoDialog";
import LabelInput from "../../labelInput";
import { labelPrompt } from "../labelUtil";

import { AnnotationsHelpers } from "../../../util/stateManager";
import actions from "../../../actions";

type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
    annotations: (state as any).annotations,
    schema: (state as any).annoMatrix?.schema,
}))
class AnnoDialogEditCategoryName extends React.PureComponent<{}, State> {
    constructor(props: {}) {
        super(props);
        this.state = {
            // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type '{... Remove this comment to see the full error message
            newCategoryText: props.metadataField,
        };
    }
    handleChangeOrSelect = (name: any) => {
        this.setState({
            newCategoryText: name,
        });
    };
    disableEditCategoryMode = () => {
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
        const { dispatch, metadataField } = this.props;
        dispatch({
            type: "annotation: disable category edit mode",
        });
        this.setState({ newCategoryText: metadataField });
    };
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
        const allCategoryNames = schema.annotations.obs.columns.map((c: any) => c.name);
        if ((allCategoryNames.indexOf(newCategoryText) > -1 &&
            newCategoryText !== metadataField) ||
            newCategoryText === "") {
            return;
        }
        this.disableEditCategoryMode();
        if (metadataField !== newCategoryText)
            dispatch(actions.annotationRenameCategoryAction(metadataField, newCategoryText));
        e.preventDefault();
    };
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
        const allCategoryNames = schema.annotations.obs.columns.map((c: any) => c.name);
        const categoryNameAlreadyExists = allCategoryNames.indexOf(name) > -1;
        const sameName = name === metadataField;
        if (categoryNameAlreadyExists && !sameName) {
            return "duplicate";
        }
        /* otherwise, no error */
        return false;
    };
    instruction = (name: any) => {
        return labelPrompt(this.editedCategoryNameError(name), "New, unique category name", ":");
    };
    allCategoryNames() {
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'schema' does not exist on type 'Readonly... Remove this comment to see the full error message
        const { schema } = this.props;
        return schema.annotations.obs.columns.map((c: any) => c.name);
    }
    render() {
        const { newCategoryText } = this.state;
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
        const { metadataField, annotations } = this.props;
        return (<>
        {/* @ts-expect-error ts-migrate(2322) FIXME: Type '{ isActive: any; inputProps: { "data-testid"... Remove this comment to see the full error message */}
        <AnnoDialog isActive={annotations.isEditingCategoryName &&
                annotations.categoryBeingEdited === metadataField} inputProps={{
                "data-testid": `${metadataField}:edit-category-name-dialog`,
            }} primaryButtonProps={{
                "data-testid": `${metadataField}:submit-category-edit`,
            // @ts-expect-error ts-migrate(2322) FIXME: Type '{ label: any; labelSuggestions: null; onChan... Remove this comment to see the full error message
            }} title="Edit category name" instruction={this.instruction(newCategoryText)} cancelTooltipContent="Close this dialog without editing this category." primaryButtonText="Edit category name" text={newCategoryText} validationError={this.editedCategoryNameError(newCategoryText)} handleSubmit={this.handleEditCategory} handleCancel={this.disableEditCategoryMode} annoInput={<LabelInput label={newCategoryText} labelSuggestions={null} onChange={this.handleChangeOrSelect} onSelect={this.handleChangeOrSelect} inputProps={{
                    "data-testid": `${metadataField}:edit-category-name-text`,
                    leftIcon: "tag",
                    intent: "none",
                    autoFocus: true,
                }} newLabelMessage="New category"/>}/>
      </>);
    }
}
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    this.state = {
    // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'props'.
    newCategoryText: (props as any).metadataField,
};
  }

  // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'handleChangeOrSelect'.
  handleChangeOrSelect = (name: any) => {
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    this.setState({
      newCategoryText: name,
    });
  };

  // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'disableEditCategoryMode'.
  disableEditCategoryMode = () => {
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "annotation: disable category edit mode",
    });
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    this.setState({ newCategoryText: metadataField });
  };

  // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'handleEditCategory'.
  handleEditCategory = (e: any) => {
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    const { dispatch, metadataField } = this.props;
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    const { newCategoryText } = this.state;

    /*
    test for uniqueness against *all* annotation names, not just the subset
    we render as categorical.
    */
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    const { schema } = this.props;
    const allCategoryNames = schema.annotations.obs.columns.map((c: any) => c.name);

    if (
      (allCategoryNames.indexOf(newCategoryText) > -1 &&
        newCategoryText !== metadataField) ||
      newCategoryText === ""
    ) {
      return;
    }

    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    this.disableEditCategoryMode();

    if (metadataField !== newCategoryText)
      dispatch(
        actions.annotationRenameCategoryAction(metadataField, newCategoryText)
      );
    e.preventDefault();
  };

  // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'editedCategoryNameError'.
  editedCategoryNameError = (name: any) => {
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
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
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    const { schema } = this.props;
    const allCategoryNames = schema.annotations.obs.columns.map((c: any) => c.name);

    const categoryNameAlreadyExists = allCategoryNames.indexOf(name) > -1;
    const sameName = name === metadataField;
    if (categoryNameAlreadyExists && !sameName) {
      return "duplicate";
    }

    /* otherwise, no error */
    return false;
  };

  // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'instruction'.
  instruction = (name: any) => {
    return labelPrompt(
      // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
      this.editedCategoryNameError(name),
      "New, unique category name",
      ":"
    );
  };

  // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'allCategoryNames'.
  allCategoryNames() {
    // @ts-expect-error ts-migrate(6133) FIXME: 'schema' is declared but its value is never read.
    const { schema } = this.props;
    return schema.annotations.obs.columns.map((c: any) => c.name);
  }

  // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'render'.
  render() {
    // @ts-expect-error ts-migrate(6133) FIXME: 'newCategoryText' is declared but its value is nev... Remove this comment to see the full error message
    const { newCategoryText } = this.state;
    // @ts-expect-error ts-migrate(6198) FIXME: All destructured elements are unused.
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
