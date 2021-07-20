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
        if (e)
            e.preventDefault();
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
        return !genesets.has(genesetName) &&
            // eslint-disable-next-line no-control-regex -- unicode 0-31 127-65535
            genesetName.match(/^\s|[\u0000-\u001F\u007F-\uFFFF]|[ ]{2,}|^$|\s$/g)
                ?.length;
    };
    render() {
        const { newGenesetName, newGenesetDescription } = this.state;
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'genesetsUI' does not exist on type 'Read... Remove this comment to see the full error message
        const { genesetsUI, parentGeneset, parentGenesetDescription, genesets, } = this.props;
        return (<>
        {/* @ts-expect-error ts-migrate(2322) FIXME: Type '{ isActive: boolean; inputProps: { "data-tes... Remove this comment to see the full error message */}
        <AnnoDialog isActive={genesetsUI.isEditingGenesetName === parentGeneset} inputProps={{
                "data-testid": `${genesetsUI.isEditingGenesetName}:rename-geneset-dialog`,
            }} primaryButtonProps={{
                "data-testid": `${genesetsUI.isEditingGenesetName}:submit-geneset`,
            }} title="Edit gene set name and description" instruction={`Rename ${genesetsUI.isEditingGenesetName}`} cancelTooltipContent="Close this dialog without renaming the gene set." primaryButtonText="Edit gene set name and description" text={newGenesetName} secondaryText={newGenesetDescription} validationError={this.validate(newGenesetName, genesets) &&
                // @ts-expect-error ts-migrate(2322) FIXME: Type '{ label: any; onChange: (e: any) => void; in... Remove this comment to see the full error message
                parentGenesetDescription === newGenesetDescription} annoInput={<LabelInput label={newGenesetName} onChange={this.handleChange} inputProps={{
                    "data-testid": "rename-geneset-modal",
                    leftIcon: "manually-entered-data",
                    intent: "none",
                    autoFocus: true,
                // @ts-expect-error ts-migrate(2322) FIXME: Type '{ label: any; onChange: (e: any) => void; in... Remove this comment to see the full error message
                }}/>} secondaryInstructions="Edit description" secondaryInput={<LabelInput label={newGenesetDescription} onChange={this.handleChangeDescription} inputProps={{ "data-testid": "change geneset description" }} intent="none" autoFocus={false}/>} handleSubmit={this.renameGeneset} handleCancel={this.disableEditGenesetNameMode}/>
      </>);
    }
}
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    this.state = {
    // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'props'.
    newGenesetName: (props as any).parentGeneset,
    // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'props'.
    newGenesetDescription: (props as any).parentGenesetDescription,
};
  }

  // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'disableEditGenesetNameMode'.
  disableEditGenesetNameMode = (e: any) => {
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    const { dispatch } = this.props;
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    this.setState({
      newGenesetName: "",
      newGenesetDescription: "",
    });
    dispatch({
      type: "geneset: disable rename geneset mode",
    });
    if (e) e.preventDefault();
  };

  // @ts-expect-error ts-migrate(2552) FIXME: Cannot find name 'renameGeneset'. Did you mean 'Re... Remove this comment to see the full error message
  renameGeneset = (e: any) => {
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    const { dispatch, genesetsUI } = this.props;
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
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

  // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'genesetNameError'.
  genesetNameError = () => {
    return false;
  };

  // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'handleChange'.
  handleChange = (e: any) => {
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    this.setState({ newGenesetName: e });
  };

  // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'handleChangeDescription'.
  handleChangeDescription = (e: any) => {
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    this.setState({ newGenesetDescription: e });
  };

  // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'validate'.
  validate = (genesetName: any, genesets: any) => {
    return !genesets.has(genesetName) &&
    // eslint-disable-next-line no-control-regex -- unicode 0-31 127-65535
    genesetName.match(/^\s|[\u0000-\u001F\u007F-\uFFFF]|[ ]{2,}|^$|\s$/g)
      ?.length;
  };

  // @ts-expect-error ts-migrate(2304) FIXME: Cannot find name 'render'.
  render() {
    // @ts-expect-error ts-migrate(6198) FIXME: All destructured elements are unused.
    const { newGenesetName, newGenesetDescription } = this.state;
    // @ts-expect-error ts-migrate(6198) FIXME: All destructured elements are unused.
    const {
      genesetsUI,
      parentGeneset,
      parentGenesetDescription,
      genesets,
    // @ts-expect-error ts-migrate(2532) FIXME: Object is possibly 'undefined'.
    } = this.props;

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
