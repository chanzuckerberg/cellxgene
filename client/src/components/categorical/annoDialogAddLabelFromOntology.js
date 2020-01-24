import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "./annoDialog";
import AnnoInputs from "./annoInputs";
import OntologySelect from "./ontologySelect";

import { labelErrorMessage, isLabelErroneous } from "./labelUtil";

@connect(state => ({
  colorAccessor: state.colors.colorAccessor,
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  universe: state.universe,
  ontology: state.ontology,
  ontologyLoading: state.ontology?.loading,
  ontologyEnabled: state.ontology?.enabled
}))
class Category extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      newLabelText: ""
    };
  }

  disableAddNewLabelFromOntologyMode = () => {
    const { dispatch } = this.props;
    dispatch({
      type: "annotation: disable add new ontology label mode"
    });
    this.setState({
      newLabelText: ""
    });
  };

  handleAddNewLabelToCategory = () => {
    const { dispatch, metadataField } = this.props;
    const { newLabelText } = this.state;

    dispatch({
      type: "annotation: add new label to category",
      metadataField,
      newLabelText,
      assignSelectedCells: false
    });
    this.setState({ newLabelText: "" });
  };

  addLabelAndAssignCells = () => {
    const { dispatch, metadataField } = this.props;
    const { newLabelText } = this.state;

    dispatch({
      type: "annotation: add new label to category",
      metadataField,
      newLabelText,
      assignSelectedCells: true
    });

    this.setState({ newLabelText: "" });
  };

  handleCreateArbitraryLabel = newLabelTextNotInOntology => {
    const { dispatch, metadataField } = this.props;

    dispatch({
      type: "annotation: add new label to category",
      metadataField,
      newLabelText: newLabelTextNotInOntology,
      assignSelectedCells: false
    });
    this.setState({ newLabelText: "" });
  };

  labelNameError = name => {
    const { metadataField, ontology, universe } = this.props;
    return isLabelErroneous(name, metadataField, ontology, universe.schema);
  };

  labelNameErrorMessage = name => {
    const { metadataField, ontology, universe } = this.props;
    return labelErrorMessage(name, metadataField, ontology, universe.schema);
  };

  /* leaky to have both of these in multiple components */
  handleChoice = e => {
    this.setState({ newLabelText: e.target });
  };

  handleTextChange = text => {
    this.setState({ newLabelText: text });
  };

  render() {
    const { newLabelText } = this.state;
    const { metadataField, annotations, ontology } = this.props;

    return (
      <>
        <AnnoDialog
          isActive={
            annotations.isAddingNewLabelFromOntology &&
            annotations.categoryAddingNewLabelFromOntology === metadataField
          }
          title="Add new label to category from existing ontology terms"
          instruction="Choose an ontology term to use as a label name:"
          cancelTooltipContent="Close this dialog without adding a label."
          primaryButtonText="Add label"
          secondaryButtonText="Add label and assign currently selected cells"
          handleSecondaryButtonSubmit={this.addLabelAndAssignCells}
          text={newLabelText}
          validationError={this.labelNameError(newLabelText)}
          errorMessage={this.labelNameErrorMessage(newLabelText)}
          handleSubmit={this.handleAddNewLabelToCategory}
          handleCancel={this.disableAddNewLabelFromOntologyMode}
          ontologySelect={
            <OntologySelect
              handleChooseOntologyTermFromDropdown={term => {
                console.log("handleChooseOntologyTermFromDropdown", term);
                this.handleChoice(term);
              }}
              categoryToDuplicate={newLabelText}
              ontology={ontology}
            />
          }
          annoInput={null}
        />
      </>
    );
  }
}

export default Category;
