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
      newGenesetName: props.originalGenesetName,
      newGenesetDescription: props.originalGenesetDescription,
      nameErrorMessage: "",
    };
  }

  disableEditGenesetNameMode = (e) => {
    const { dispatch } = this.props;
    const { originalGenesetDescription, originalGenesetName } = this.props;
    this.setState({
      newGenesetName: originalGenesetName,
      newGenesetDescription: originalGenesetDescription,
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

  genesetNameError = () => false;

  handleChange = (e) => {
    this.setState({ newGenesetName: e });
  };

  handleChangeDescription = (e) => {
    this.setState({ newGenesetDescription: e });
  };

  validate = (
    originalGenesetName,
    newGenesetName,
    originalGenesetDescription,
    newGenesetDescription,
    genesets
  ) => {
    if (
      originalGenesetName !== newGenesetName &&
      genesets.has(newGenesetName)
    ) {
      this.setState({
        nameErrorMessage: `There is already a gene set with that name`,
      });
      return true;
    }
    if (
      // eslint-disable-next-line no-control-regex -- unicode 0-31 127-65535
      newGenesetName.match(/^[\u0000-\u001F\u007F-\uFFFF]|[ ]{2,}/g)?.length
    ) {
      this.setState({
        nameErrorMessage:
          "Gene set names can only contain alphanumeric characters and the following special characters: ! ” # $ % ’ ( ) * + , - . / : ; < = > ? @  ] ^ _ ` | ~",
      });
      return true;
    }
    this.setState({
      nameErrorMessage: "",
    });
    if (
      originalGenesetName === newGenesetName &&
      originalGenesetDescription === newGenesetDescription
    )
      return true;
    return false;
  };

  render() {
    const { newGenesetName, newGenesetDescription, nameErrorMessage } =
      this.state;
    const {
      genesetsUI,
      originalGenesetName,
      originalGenesetDescription,
      genesets,
    } = this.props;

    return (
      <AnnoDialog
          isActive={genesetsUI.isEditingGenesetName === originalGenesetName}
          inputProps={{
            "data-testid": `${genesetsUI.isEditingGenesetName}:rename-geneset-dialog`,
          }}
          primaryButtonProps={{
            "data-testid": `${genesetsUI.isEditingGenesetName}:submit-geneset`,
          }}
          title="Edit gene set name and description"
          instruction={
            <>
              Rename <b>{genesetsUI.isEditingGenesetName}</b>
            </>
          }
          cancelTooltipContent="Close this dialog without renaming the gene set."
          primaryButtonText="Edit gene set name and description"
          text={newGenesetName}
          secondaryText={newGenesetDescription}
          validationError={this.validate(
            originalGenesetName,
            newGenesetName,
            originalGenesetDescription,
            newGenesetDescription,
            genesets
          )}
          errorMessage={nameErrorMessage}
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
    );
  }
}

export default RenameGeneset;
