import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../annoDialog";
import LabelInput from "../../labelInput";

type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  annotations: (state as any).annotations,
  schema: (state as any).annoMatrix?.schema,
  obsCrossfilter: (state as any).obsCrossfilter,
  genesetsUI: (state as any).genesetsUI,
  genesets: (state as any).genesets.genesets,
}))
class RenameGeneset extends React.PureComponent<{}, State> {
  constructor(props: {}) {
    super(props);
    this.state = {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'parentGeneset' does not exist on type '{... Remove this comment to see the full error message
      newGenesetName: props.parentGeneset,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'parentGenesetDescription' does not exist... Remove this comment to see the full error message
      newGenesetDescription: props.parentGenesetDescription,
    };
  }

  disableEditGenesetNameMode = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
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

  renameGeneset = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
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

  handleChange = (e: any) => {
    this.setState({ newGenesetName: e });
  };

  handleChangeDescription = (e: any) => {
    this.setState({ newGenesetDescription: e });
  };

  validate = (genesetName: any, genesets: any) => {
    return (
      !genesets.has(genesetName) &&
      // eslint-disable-next-line no-control-regex -- unicode 0-31 127-65535
      genesetName.match(/^\s|[\u0000-\u001F\u007F-\uFFFF]|[ ]{2,}|^$|\s$/g)
        ?.length
    );
  };

  render() {
    const { newGenesetName, newGenesetDescription } = this.state;
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'genesetsUI' does not exist on type 'Read... Remove this comment to see the full error message
      genesetsUI,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'parentGeneset' does not exist on type 'Read... Remove this comment to see the full error message
      parentGeneset,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'parentGenesetDescription' does not exist on type 'Read... Remove this comment to see the full error message
      parentGenesetDescription,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'genesets' does not exist on type 'Read... Remove this comment to see the full error message
      genesets,
    } = this.props;
    return (
      <>
        <AnnoDialog
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ isActive: boolean; inputProps: { "data-tes... Remove this comment to see the full error message
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
              // @ts-expect-error ts-migrate(2322) FIXME: Type '{ label: any; onChange: (e: any) => void; in... Remove this comment to see the full error message
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
              // @ts-expect-error ts-migrate(2322) FIXME: Type '{ label: any; onChange: (e: any) => void; in... Remove this comment to see the full error message
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
