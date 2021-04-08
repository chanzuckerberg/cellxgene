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
      newGenesetName: props.parentGeneset,
      newGenesetDescription: props.parentGenesetDescription,
    };
  }

  disableEditGenesetNameMode = (e) => {
    const { dispatch } = this.props;
    this.setState({
      newGenesetName: "",
      newGenesetDescription: "",
    });
    dispatch({
      type: "geneset: disable rename geneset mode",
    });
    if (e) e.preventDefault();
  };

  renameGeneset = (e) => {
    const { dispatch, genesetsUI } = this.props;
    const { newGenesetName, newGenesetDescription } = this.state;

    dispatch({
      type: "geneset: update",
      genesetName: genesetsUI.isEditingGenesetName,
      update: {
        genesetName: newGenesetName,
        genesetDescription: newGenesetDescription,
      },
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
    this.setState({ newGenesetName: e });
  };

  handleChangeDescription = (e) => {
    this.setState({ newGenesetDescription: e });
  };

  render() {
    const { newGenesetName, newGenesetDescription } = this.state;
    const { genesetsUI, parentGeneset, parentGenesetDescription } = this.props;

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
          title="Edit gene set and description"
          instruction={`Rename ${genesetsUI.isEditingGenesetName}`}
          cancelTooltipContent="Close this dialog without renaming the gene set."
          primaryButtonText="Edit gene set and description"
          text={newGenesetName}
          secondaryText={newGenesetDescription}
          validationError={
            genesetsUI.isEditingGenesetName === newGenesetName &&
            parentGenesetDescription === newGenesetDescription
          }
          annoInput={
            <LabelInput
              label={newGenesetName}
              onChange={this.handleChange}
              inputProps={{
                "data-testid": "rename-geneset-modal",
                leftIcon: "manually-entered-data",
                intent: "none",
                autoFocus: true,
              }}
            />
          }
          secondaryInstructions="Edit description"
          secondaryInput={
            <LabelInput
              label={newGenesetDescription}
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
