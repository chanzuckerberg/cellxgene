import _ from "lodash";
import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../annoDialog";
import LabelInput from "../../labelInput";

@connect((state) => ({
  annotations: state.annotations,
  schema: state.annoMatrix?.schema,
  ontology: state.ontology,
  obsCrossfilter: state.obsCrossfilter,
  genesets: state.genesets.genesets,
  genesetsUI: state.genesetsUI,
}))
class CreateGenesetDialogue extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      genesetName: "",
      genesToPopulateGeneset: "",
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
      genesetName,
      genesetDescription: "",
    });
    if (genesToPopulateGeneset) {
      const genesTmpHardcodedFormat = [];

      const genesArrayFromString = _.pull(
        _.uniq(genesToPopulateGeneset.split(/[ ,]+/)),
        ""
      );

      genesArrayFromString.forEach((_gene) => {
        genesTmpHardcodedFormat.push({
          geneSymbol: _gene,
        });
      });

      dispatch({
        type: "geneset: add genes",
        genesetName,
        genes: genesTmpHardcodedFormat,
      });
    }
    dispatch({
      type: "geneset: disable create geneset mode",
    });
    this.setState({
      genesetName: "",
      genesToPopulateGeneset: "",
    });
    e.preventDefault();
  };

  genesetNameError = () => {
    return false;
  };

  handleChange = (e) => {
    this.setState({ genesetName: e });
  };

  handleGenesetInputChange = (e) => {
    this.setState({ genesToPopulateGeneset: e });
  };

  instruction = (genesetName, genesets) => {
    return genesets.has(genesetName)
      ? "Geneset name must be unique."
      : "New, unique geneset name";
  };

  validate = (genesetName, genesets) => {
    return genesets.has(genesetName);
  };

  render() {
    const { genesetName } = this.state;
    const { metadataField, genesetsUI, genesets } = this.props;

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
          instruction={this.instruction(genesetName, genesets)}
          cancelTooltipContent="Close this dialog without creating a new gene set."
          primaryButtonText="Create gene set"
          text={genesetName}
          validationError={this.validate(genesetName, genesets)}
          annoInput={
            <LabelInput
              onChange={this.handleChange}
              inputProps={{
                "data-testid": "create-geneset-modal",
                leftIcon: "manually-entered-data",
                intent: "none",
                autoFocus: true,
              }}
              newLabelMessage="Create gene set"
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
