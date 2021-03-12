import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../categorical/annoDialog";
import LabelInput from "../../categorical/labelInput";

@connect((state) => ({
  annotations: state.annotations,
  schema: state.annoMatrix?.schema,
  ontology: state.ontology,
  obsCrossfilter: state.obsCrossfilter,
  genesetsUI: state.genesetsUI,
}))
class RenameGeneset extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      genesetName: "",
    };
  }

  disableEditGenesetNameMode = (e) => {
    const { dispatch } = this.props;
    this.setState({
      genesetName: "",
    });
    dispatch({
      type: "geneset: disable rename geneset mode",
    });
    if (e) e.preventDefault();
  };

  renameGeneset = (e) => {
    const { dispatch, genesetsUI } = this.props;
    const { genesetName } = this.state;

    dispatch({
      type: "geneset: rename",
      name: genesetsUI.isEditingGenesetName,
      newName: genesetName,
    });
    dispatch({
      type: "geneset: disable rename geneset mode",
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

  instruction = () => {
    return "New, unique geneset name";
    /* todo genesets */
    // return genesetPrompt(this.genesetNameError(geneset), "New, unique geneset", ":");
  };

  render() {
    const { genesetName } = this.state;
    const { genesetsUI } = this.props;

    return (
      <>
        <AnnoDialog
          isActive={genesetsUI.isEditingGenesetName}
          inputProps={{
            "data-testid": `${genesetsUI.isEditingGenesetName}:rename-geneset-dialog`,
          }}
          primaryButtonProps={{
            "data-testid": `${genesetsUI.isEditingGenesetName}:submit-geneset`,
          }}
          title="Rename gene set"
          instruction={`Rename ${genesetsUI.isEditingGenesetName}`} /* todo genesets this.instruction(genesetName) */
          cancelTooltipContent="Close this dialog without renaming the gene set."
          primaryButtonText="Rename gene set"
          text={genesetName}
          validationError={false}
          annoInput={
            <LabelInput
              onChange={this.handleChange}
              inputProps={{
                "data-testid": "rename-geneset-modal",
                leftIcon: "manually-entered-data",
                intent: "none",
                autoFocus: true,
              }}
              newLabelMessage="Rename gene set"
            />
          }
          handleSubmit={this.renameGeneset}
          handleCancel={this.disableEditGenesetNameMode}
        />
      </>
    );
  }
}

export default RenameGeneset;
