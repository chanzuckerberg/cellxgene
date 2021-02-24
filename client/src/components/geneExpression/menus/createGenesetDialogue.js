import _ from "lodash";
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
      genesToPopulateGeneset: "",
    });
    dispatch({
      type: "geneset: disable create geneset mode",
    });
    if (e) e.preventDefault();
  };

  createGeneset = (e) => {
    const { dispatch } = this.props;
    const { genesetName, genesToPopulateGeneset } = this.state;

    dispatch({
      type: "geneset: create",
      name: genesetName,
      genes: genesToPopulateGeneset
        ? _.pull(_.uniq(genesToPopulateGeneset.split(/[ ,]+/)), "")
        : null,
    });
    dispatch({
      type: "geneset: disable create geneset mode",
    });
    e.preventDefault();
  };

  genesetNameError = () => {
    /* todo genesets validation */
    return false;
  };

  handleChange = (e) => {
    this.setState({ genesetName: e });
  };

  handleGenesetInputChange = (e) => {
    this.setState({ genesToPopulateGeneset: e });
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
          inputProps={{
            "data-testid": `${metadataField}:create-geneset-dialog`,
          }}
          primaryButtonProps={{
            "data-testid": `${metadataField}:submit-geneset`,
          }}
          title="Create gene set"
          instruction="Type a new name for your gene set" /* todo genesets this.instruction(genesetName) */
          cancelTooltipContent="Close this dialog without creating the gene set."
          primaryButtonText="Create gene set same"
          text={genesetName}
          validationError={false}
          annoInput={
            <LabelInput
              onChange={this.handleChange}
              inputProps={{
                "data-testid": "create-geneset-modal",
                leftIcon: "manually-entered-data",
                intent: "none",
                autoFocus: true,
              }}
              newLabelMessage="Create gene set name"
            />
          }
          secondaryInstructions="Optionally add a list of comma separated genes to populate the gene set"
          secondaryInput={
            <LabelInput
              onChange={this.handleGenesetInputChange}
              inputProps={{
                "data-testid": "add-genes",
                intent: "none",
                autoFocus: false,
              }}
              newLabelMessage="populate geneset with genes"
            />
          }
          handleSubmit={this.createGeneset}
          handleCancel={this.disableCreateGenesetMode}
        />
      </>
    );
  }
}

export default CreateGenesetDialogue;
