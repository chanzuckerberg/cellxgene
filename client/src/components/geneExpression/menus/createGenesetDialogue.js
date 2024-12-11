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
      type: "sampleset: disable create sampleset mode",
    });
    if (e) e.preventDefault();
  };

  createGeneset = (e) => {
    const { dispatch } = this.props;
    const { genesetName, genesToPopulateGeneset, genesetDescription } =
      this.state;

    dispatch({
      type: "sampleset: create",
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
      type: "sampleset: disable create sampleset mode",
    });
    this.setState({
      genesetName: "",
      genesToPopulateGeneset: "",
    });
    e.preventDefault();
  };

  genesetNameError = () => false;

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

  instruction = (genesetName, genesets) =>
    genesets.has(genesetName)
      ? "Sample set name must be unique."
      : "New, unique sample set name";

  validate = (genesetName, genesets) => {
    if (genesets.has(genesetName)) {
      this.setState({
        nameErrorMessage: "There is already a sample set with that name",
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
          "Sample set names can only contain alphanumeric characters and the following special characters: ! ” # $ % ’ ( ) * + , - . / : ; < = > ? @  ] ^ _ ` | ~",
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
    const { genesetsUI, genesets } = this.props;

    return (
      <>
        <Dialog
          icon="tag"
          title="Create sample set"
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
                    "data-testid": "create-geneset-input",
                    leftIcon: "manually-entered-data",
                    intent: "none",
                    autoFocus: true,
                  }}
                  newLabelMessage="Create sample set"
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
                  sample set
                </p>
                <LabelInput
                  onChange={this.handleDescriptionInputChange}
                  inputProps={{
                    "data-testid": "add-geneset-description",
                    intent: "none",
                    autoFocus: false,
                  }}
                  newLabelMessage="Add sampleset description"
                />

                <p style={{ marginTop: 20 }}>
                  Optionally add a list of comma separated{" "}
                  <span style={{ fontWeight: 700 }}>samples</span> to populate
                  the sample set
                </p>
                <LabelInput
                  onChange={this.handleGenesetInputChange}
                  inputProps={{
                    "data-testid": "add-genes",
                    intent: "none",
                    autoFocus: false,
                  }}
                  newLabelMessage="populate sampleset with samples"
                />
              </div>
            </div>
            <div className={Classes.DIALOG_FOOTER}>
              <div className={Classes.DIALOG_FOOTER_ACTIONS}>
                <Tooltip2 content="Close this dialog without creating a new sample set.">
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
                  Create sample set
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
