import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../annoDialog";
import LabelInput from "../../labelInput";

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
      genesetName: props.parentGeneset,
      genesetDescription: props.parentGenesetDescription,
    };
  }

  disableEditGenesetNameMode = (e) => {
    const { dispatch } = this.props;
    this.setState({
      genesetName: "",
      genesetDescription: "",
    });
    dispatch({
      type: "geneset: disable rename geneset mode",
    });
    if (e) e.preventDefault();
  };

  renameGeneset = (e) => {
    const { dispatch, genesetsUI } = this.props;
    const { genesetName, genesetDescription } = this.state;

    dispatch({
      type: "geneset: update",
      genesetName: genesetsUI.isEditingGenesetName,
      update: { genesetName, genesetDescription },
    });
    dispatch({
      type: "geneset: disable rename geneset mode",
    });
    e.preventDefault();
  };

  genesetNameError = () => {
    return false;
  };

  handleChange = (e) => {
    this.setState({ genesetName: e });
  };

  handleChangeDescription = (e) => {
    this.setState({ genesetDescription: e });
  };

  render() {
    const { genesetName, genesetDescription } = this.state;
    const { genesetsUI, parentGeneset } = this.props;

    return (
      <>
        <AnnoDialog
          isActive={genesetsUI.isEditingGenesetName === parentGeneset}
          inputProps={{
            "data-testid": `${genesetsUI.isEditingGenesetName}:rename-geneset-dialog`,
          }}
          primaryButtonProps={{
            "data-testid": `${genesetsUI.isEditingGenesetName}:submit-geneset`,
          }}
          title="Rename gene set"
          instruction={`Rename ${genesetsUI.isEditingGenesetName}`}
          cancelTooltipContent="Close this dialog without renaming the gene set."
          primaryButtonText="Rename gene set"
          text={genesetName}
          secondaryText={genesetDescription}
          validationError={genesetsUI.isEditingGenesetName === genesetName}
          annoInput={
            <LabelInput
              label={genesetName}
              onChange={this.handleChange}
              inputProps={{
                "data-testid": "rename-geneset-modal",
                leftIcon: "manually-entered-data",
                intent: "none",
                autoFocus: true,
              }}
            />
          }
          secondaryInstructions="Geneset Description"
          secondaryInput={
            <LabelInput
              label={genesetDescription}
              onChange={this.handleChangeDescription}
              inputProps={{ "data-testid": "change geneset description" }}
              intent="none"
              autoFocus={false}
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
