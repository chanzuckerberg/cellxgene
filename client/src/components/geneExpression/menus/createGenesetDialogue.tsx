import pull from "lodash.pull";
import uniq from "lodash.uniq";
import React from "react";
import { connect } from "react-redux";
import { Button, Dialog, Classes, Colors } from "@blueprintjs/core";
import { Tooltip2 } from "@blueprintjs/popover2";
import LabelInput from "../../labelInput";
import actions from "../../../actions";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  annotations: (state as any).annotations,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  schema: (state as any).annoMatrix?.schema,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  obsCrossfilter: (state as any).obsCrossfilter,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  genesets: (state as any).genesets.genesets,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  genesetsUI: (state as any).genesetsUI,
}))
// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
class CreateGenesetDialogue extends React.PureComponent<{}, State> {
  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);
    this.state = {
      genesetName: "",
      genesToPopulateGeneset: "",
      genesetDescription: "",
    };
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  disableCreateGenesetMode = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    this.setState({
      genesetName: "",
      genesToPopulateGeneset: "",
      genesetDescription: "",
      nameErrorMessage: "",
    });
    dispatch({
      type: "geneset: disable create geneset mode",
    });
    if (e) e.preventDefault();
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  createGeneset = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    const {
      genesetName,
      genesToPopulateGeneset,
      genesetDescription,
    } = this.state;
    dispatch({
      type: "geneset: create",
      genesetName: genesetName.trim(),
      genesetDescription,
    });
    if (genesToPopulateGeneset) {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      const genesTmpHardcodedFormat: any = [];
      const genesArrayFromString = pull(
        uniq(genesToPopulateGeneset.split(/[ ,]+/)),
        ""
      );
      genesArrayFromString.forEach((_gene) => {
        genesTmpHardcodedFormat.push({
          geneSymbol: _gene,
        });
      });
      dispatch(actions.genesetAddGenes(genesetName, genesTmpHardcodedFormat));
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  genesetNameError = () => {
    return false;
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleChange = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'genesets' does not exist on type 'Readon... Remove this comment to see the full error message
    const { genesets } = this.props;
    this.setState({ genesetName: e });
    this.validate(e, genesets);
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleGenesetInputChange = (e: any) => {
    this.setState({ genesToPopulateGeneset: e });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleDescriptionInputChange = (e: any) => {
    this.setState({ genesetDescription: e });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  instruction = (genesetName: any, genesets: any) => {
    return genesets.has(genesetName)
      ? "Gene set name must be unique."
      : "New, unique gene set name";
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  validate = (genesetName: any, genesets: any) => {
    if (genesets.has(genesetName)) {
      this.setState({
        nameErrorMessage: "There is already a geneset with that name",
      });
      return false;
    }
    if (
      genesetName.length > 1 &&
      // eslint-disable-next-line no-control-regex -- unicode 0-31 127-65535
      genesetName.match(/^[\u0000-\u001F\u007F-\uFFFF]|[ ]{2,}/g)?.length
    ) {
      this.setState({
        nameErrorMessage:
          "Gene set names can only contain alphanumeric characters and the following special characters: ! ” # $ % ’ ( ) * + , - . / : ; < = > ? @  ] ^ _ ` | ~",
      });
      return false;
    }
    this.setState({
      nameErrorMessage: "",
    });
    return true;
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const { genesetName, nameErrorMessage } = this.state;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'genesetsUI' does not exist on type 'Read... Remove this comment to see the full error message
    const { genesetsUI, genesets } = this.props;
    return (
      <>
        <Dialog
          icon="tag"
          title="Create gene set"
          isOpen={genesetsUI.createGenesetModeActive}
          onClose={this.disableCreateGenesetMode}
        >
          <form
            onSubmit={(e) => {
              e.preventDefault();
            }}
          >
            <div className={Classes.DIALOG_BODY}>
              <div style={{ marginBottom: 20 }}>
                <p>{this.instruction(genesetName, genesets)}</p>
                <LabelInput
                  // @ts-expect-error ts-migrate(2322) FIXME: Type '{ onChange: (e: any) => void; inputProps: { ... Remove this comment to see the full error message
                  onChange={this.handleChange}
                  inputProps={{
                    "data-testid": "create-geneset-input",
                    leftIcon: "manually-entered-data",
                    intent: "none",
                    autoFocus: true,
                  }}
                  newLabelMessage="Create gene set"
                />
                <p
                  style={{
                    marginTop: 7,
                    visibility: nameErrorMessage !== "" ? "visible" : "hidden",
                    color: Colors.ORANGE3,
                  }}
                >
                  {nameErrorMessage}
                </p>
                <p style={{ marginTop: 20 }}>
                  Optionally add a{" "}
                  <span style={{ fontWeight: 700 }}>description</span> for this
                  gene set
                </p>
                <LabelInput
                  // @ts-expect-error ts-migrate(2322) FIXME: Type '{ onChange: (e: any) => void; inputProps: { ... Remove this comment to see the full error message
                  onChange={this.handleDescriptionInputChange}
                  inputProps={{
                    "data-testid": "add-geneset-description",
                    intent: "none",
                    autoFocus: false,
                  }}
                  newLabelMessage="Add geneset description"
                />

                <p style={{ marginTop: 20 }}>
                  Optionally add a list of comma separated{" "}
                  <span style={{ fontWeight: 700 }}>genes</span> to populate the
                  gene set
                </p>
                <LabelInput
                  // @ts-expect-error ts-migrate(2322) FIXME: Type '{ onChange: (e: any) => void; inputProps: { ... Remove this comment to see the full error message
                  onChange={this.handleGenesetInputChange}
                  inputProps={{
                    "data-testid": "add-genes",
                    intent: "none",
                    autoFocus: false,
                  }}
                  newLabelMessage="populate geneset with genes"
                />
              </div>
            </div>
            <div className={Classes.DIALOG_FOOTER}>
              <div className={Classes.DIALOG_FOOTER_ACTIONS}>
                <Tooltip2 content="Close this dialog without creating a new gene set.">
                  <Button onClick={this.disableCreateGenesetMode}>
                    Cancel
                  </Button>
                </Tooltip2>
                <Button
                  data-testid="submit-geneset"
                  onClick={this.createGeneset}
                  disabled={nameErrorMessage !== ""}
                  intent="primary"
                  type="submit"
                >
                  Create gene set
                </Button>
              </div>
            </div>
          </form>
        </Dialog>
      </>
    );
  }
}

export default CreateGenesetDialogue;
