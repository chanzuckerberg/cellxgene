import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../annoDialog";
import LabelInput from "../../labelInput";
import { labelPrompt, isLabelErroneous } from "../labelUtil";
import actions from "../../../actions";

type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  annotations: (state as any).annotations,
  schema: (state as any).annoMatrix?.schema,
  obsCrossfilter: (state as any).obsCrossfilter,
}))
class Category extends React.PureComponent<{}, State> {
  constructor(props: {}) {
    super(props);
    this.state = {
      newLabelText: "",
    };
  }

  disableAddNewLabelMode = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    this.setState({
      newLabelText: "",
    });
    dispatch({
      type: "annotation: disable add new label mode",
    });
    if (e) e.preventDefault();
  };

  handleAddNewLabelToCategory = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField } = this.props;
    const { newLabelText } = this.state;
    // @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 0.
    this.disableAddNewLabelMode();
    dispatch(
      actions.annotationCreateLabelInCategory(
        metadataField,
        newLabelText,
        false
      )
    );
    e.preventDefault();
  };

  addLabelAndAssignCells = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField } = this.props;
    const { newLabelText } = this.state;
    // @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 0.
    this.disableAddNewLabelMode();
    dispatch(
      actions.annotationCreateLabelInCategory(metadataField, newLabelText, true)
    );
    e.preventDefault();
  };

  labelNameError = (name: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
    const { metadataField, schema } = this.props;
    return isLabelErroneous(name, metadataField, schema);
  };

  instruction = (label: any) => {
    return labelPrompt(this.labelNameError(label), "New, unique label", ":");
  };

  handleChangeOrSelect = (label: any) => {
    this.setState({ newLabelText: label });
  };

  render() {
    const { newLabelText } = this.state;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
    const { metadataField, annotations, obsCrossfilter } = this.props;
    return (
      <>
        <AnnoDialog
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ isActive: any; inputProps: { "data-testid"... Remove this comment to see the full error message
          isActive={
            annotations.isAddingNewLabel &&
            annotations.categoryAddingNewLabel === metadataField
          }
          inputProps={{ "data-testid": `${metadataField}:create-label-dialog` }}
          primaryButtonProps={{
            "data-testid": `${metadataField}:submit-label`,
          }}
          title="Add new label to category"
          instruction={this.instruction(newLabelText)}
          cancelTooltipContent="Close this dialog without adding a label."
          primaryButtonText="Add label"
          secondaryButtonText={`Add label & assign ${obsCrossfilter.countSelected()} selected cells`}
          handleSecondaryButtonSubmit={this.addLabelAndAssignCells}
          text={newLabelText}
          validationError={this.labelNameError(newLabelText)}
          handleSubmit={this.handleAddNewLabelToCategory}
          handleCancel={this.disableAddNewLabelMode}
          annoInput={
            <LabelInput
              // @ts-expect-error ts-migrate(2322) FIXME: Type '{ labelSuggestions: null; onChange: (label: ... Remove this comment to see the full error message
              labelSuggestions={null}
              onChange={this.handleChangeOrSelect}
              onSelect={this.handleChangeOrSelect}
              inputProps={{
                "data-testid": `${metadataField}:new-label-name`,
                leftIcon: "tag",
                intent: "none",
                autoFocus: true,
              }}
            />
          }
        />
      </>
    );
  }
}

export default Category;
