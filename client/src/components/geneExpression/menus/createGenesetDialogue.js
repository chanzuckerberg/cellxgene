import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../categorical/annoDialog";
import LabelInput from "../../categorical/labelInput";
// import LabelInput from "../labelInput";
// import { genesetPrompt, isGenesetErroneous } from "../genesetUtil";

@connect((state) => ({
  annotations: state.annotations,
  schema: state.annoMatrix?.schema,
  ontology: state.ontology,
  obsCrossfilter: state.obsCrossfilter,
  genesetsUI: state.genesetsUI,
}))
class CreateGenesetDialogue extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      genesetName: "",
    };
  }

  disableCreateGenesetMode = (e) => {
    const { dispatch } = this.props;
    this.setState({
      genesetName: "",
    });
    dispatch({
      type: "geneset: disable create geneset mode",
    });
    if (e) e.preventDefault();
  };

  createGeneset = (e) => {
    const { dispatch } = this.props;
    const { genesetName } = this.state;

    dispatch({ type: "geneset: create", genesetName });
    e.preventDefault();
  };

  genesetNameError = () => {
    /* todo genesets validation */
    return false;
  };

  handleChange = () => {
    return () => {};
  };

  handleSelect = () => {
    return () => {};
  };

  instruction = () => {
    return "New, unique geneset";
    /* todo genesets */
    // return genesetPrompt(this.genesetNameError(geneset), "New, unique geneset", ":");
  };

  render() {
    const { genesetName } = this.state;
    const { metadataField, genesetsUI } = this.props;

    return (
      <>
        <AnnoDialog
          isActive={genesetsUI.createGenesetModeActive}
          inputProps={{ "data-testid": `${metadataField}:create-label-dialog` }}
          primaryButtonProps={{
            "data-testid": `${metadataField}:submit-label`,
          }}
          title="Create new gene set"
          instruction="Give your geneset a name" /* todo genesets this.instruction(genesetName) */
          cancelTooltipContent="Close this dialog without creating a gene set."
          primaryButtonText="Create geneset"
          handleSecondaryButtonSubmit={this.addLabelAndAssignCells}
          text={genesetName}
          validationError={() => {
            /* todo genesets this.labelNameError(genesetName) */
            return false;
          }}
          annoInput={
            <LabelInput
              onChange={this.handleChange}
              onSelect={this.handleSelect}
              inputProps={{
                "data-testid": "add-genes",
                leftIcon: "manually-entered-data",
                intent: "none",
                autoFocus: true,
              }}
              newLabelMessage="New category"
            />
          }
          handleSubmit={this.handleAddNewLabelToCategory}
          handleCancel={this.disableCreateGenesetMode}
        />
      </>
    );
  }
}

export default CreateGenesetDialogue;
