import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../annoDialog";
import LabelInput from "../labelInput";
import { labelPrompt, isLabelErroneous } from "../labelUtil";

@connect((state) => ({
  colorAccessor: state.colors.colorAccessor,
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  universe: state.universe,
  ontology: state.ontology,
  crossfilter: state.crossfilter,
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
    dispatch({
      type: "annotation: add new label to category",
      metadataField,
      newLabelText,
      assignSelectedCells: false,
    });
    e.preventDefault();
  };

  addLabelAndAssignCells = (e) => {
    const { dispatch, metadataField } = this.props;
    const { newLabelText } = this.state;

    this.disableAddNewLabelMode();
    dispatch({
      type: "annotation: add new label to category",
      metadataField,
      newLabelText,
      assignSelectedCells: true,
    });
    e.preventDefault();
  };

  labelNameError = (name) => {
    const { metadataField, ontology, universe } = this.props;
    return isLabelErroneous(name, metadataField, ontology, universe.schema);
  };

  instruction = (label) => {
    return labelPrompt(this.labelNameError(label), "New, unique label", ":");
  };

  handleChangeOrSelect = (label) => {
    this.setState({ newLabelText: label });
  };

  render() {
    const { newLabelText } = this.state;
    const { metadataField, annotations, ontology, crossfilter } = this.props;
    const ontologyEnabled = ontology?.enabled ?? false;

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
          secondaryButtonText={`Add label and assign ${crossfilter.countSelected()} currently selected cells`}
          handleSecondaryButtonSubmit={this.addLabelAndAssignCells}
          text={newLabelText}
          validationError={this.labelNameError(newLabelText)}
          handleSubmit={this.handleAddNewLabelToCategory}
          handleCancel={this.disableAddNewLabelMode}
          annoInput={
            <LabelInput
              labelSuggestions={ontologyEnabled ? ontology.terms : null}
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
