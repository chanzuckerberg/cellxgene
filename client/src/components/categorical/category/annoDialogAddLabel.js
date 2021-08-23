import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../annoDialog";
import LabelInput from "../../labelInput";
import { labelPrompt, isLabelErroneous } from "../labelUtil";
import actions from "../../../actions";

@connect((state) => ({
  annotations: state.annotations,
  schema: state.annoMatrix?.schema,
  obsCrossfilter: state.obsCrossfilter,
}))
class Category extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      newLabelText: "",
    };
  }

  disableAddNewLabelMode = (e) => {
    const { dispatch } = this.props;
    this.setState({
      newLabelText: "",
    });
    dispatch({
      type: "annotation: disable add new label mode",
    });
    if (e) e.preventDefault();
  };

  handleAddNewLabelToCategory = (e) => {
    const { dispatch, metadataField } = this.props;
    const { newLabelText } = this.state;

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

  addLabelAndAssignCells = (e) => {
    const { dispatch, metadataField } = this.props;
    const { newLabelText } = this.state;

    this.disableAddNewLabelMode();
    dispatch(
      actions.annotationCreateLabelInCategory(metadataField, newLabelText, true)
    );
    e.preventDefault();
  };

  labelNameError = (name) => {
    const { metadataField, schema } = this.props;
    return isLabelErroneous(name, metadataField, schema);
  };

  instruction = (label) => labelPrompt(this.labelNameError(label), "New, unique label", ":");

  handleChangeOrSelect = (label) => {
    this.setState({ newLabelText: label });
  };

  render() {
    const { newLabelText } = this.state;
    const { metadataField, annotations, obsCrossfilter } = this.props;

    return (
      <>
        <AnnoDialog
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
