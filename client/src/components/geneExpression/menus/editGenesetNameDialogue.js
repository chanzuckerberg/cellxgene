import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../annoDialog";
import LabelInput from "../../labelInput";

@connect((state) => ({
  annotations: state.annotations,
  schema: state.annoMatrix?.schema,
  obsCrossfilter: state.obsCrossfilter,
  genesetsUI: state.genesetsUI,
  genesets: state.genesets.genesets,
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

  validate = (genesetName, genesets) => {
    return (
      !genesets.has(genesetName) &&
      // eslint-disable-next-line no-control-regex -- unicode 0-31 127-65535
      genesetName.match(/^\s|[\u0000-\u001F\u007F-\uFFFF]|[ ]{2,}|^$|\s$/g)
        ?.length
    );
  };

  render() {
    const { newGenesetName, newGenesetDescription } = this.state;
    const { genesetsUI, parentGeneset, parentGenesetDescription, genesets } =
      this.props;

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
          title="Edit gene set name and description"
          instruction={`Rename ${genesetsUI.isEditingGenesetName}`}
          cancelTooltipContent="Close this dialog without renaming the gene set."
          primaryButtonText="Edit gene set name and description"
          text={newGenesetName}
          secondaryText={newGenesetDescription}
          validationError={
            this.validate(newGenesetName, genesets) &&
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
