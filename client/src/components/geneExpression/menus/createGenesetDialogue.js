import pull from "lodash.pull";
import uniq from "lodash.uniq";
import React from "react";
import { connect } from "react-redux";
import { Button, Dialog, Classes, Colors } from "@blueprintjs/core";
import { Tooltip2 } from "@blueprintjs/popover2";
import LabelInput from "../../labelInput";
import actions from "../../../actions";

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
      genesetDescription: "",
    };
  }

  disableCreateGenesetMode = (e) => {
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

  createGeneset = (e) => {
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
      const genesTmpHardcodedFormat = [];

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

  genesetNameError = () => {
    return false;
  };

  handleChange = (e) => {
    const { genesets } = this.props;
    this.setState({ genesetName: e });
    this.validate(e, genesets);
  };

  handleGenesetInputChange = (e) => {
    this.setState({ genesToPopulateGeneset: e });
  };

  handleDescriptionInputChange = (e) => {
    this.setState({ genesetDescription: e });
  };

  instruction = (genesetName, genesets) => {
    return genesets.has(genesetName)
      ? "Gene set name must be unique."
      : "New, unique gene set name";
  };

  validate = (genesetName, genesets) => {
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

  render() {
    const { genesetName, nameErrorMessage } = this.state;
    const { metadataField, genesetsUI, genesets } = this.props;

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
                  onChange={this.handleChange}
                  inputProps={{
                    "data-testid": "create-geneset-modal",
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
                  data-testid={`${metadataField}:submit-geneset`}
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
