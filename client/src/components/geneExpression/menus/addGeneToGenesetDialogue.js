import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../annoDialog";
import LabelInput from "../../labelInput";
import parseBulkGeneString from "../../../util/parseBulkGeneString";
import actions from "../../../actions";

@connect((state) => ({
  genesetsUI: state.genesetsUI,
}))
class AddGeneToGenesetDialogue extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      genesToAdd: "",
    };
  }

  disableAddGeneMode = () => {
    const { dispatch } = this.props;
    dispatch({
      type: "sampleset: disable add new samples mode",
    });
  };

  handleAddGeneToGeneSet = (e) => {
    const { geneset, dispatch } = this.props;
    const { genesToAdd } = this.state;

    const genesTmpHardcodedFormat = [];

    const genesArrayFromString = parseBulkGeneString(genesToAdd);

    genesArrayFromString.forEach((_gene) => {
      genesTmpHardcodedFormat.push({
        geneSymbol: _gene,
      });
    });

    dispatch(actions.genesetAddGenes(geneset, genesTmpHardcodedFormat));
    dispatch({
      type: "sampleset: disable add new samples mode",
    });
    if (e) e.preventDefault();
  };

  handleChange = (e) => {
    this.setState({ genesToAdd: e });
  };

  render() {
    const { geneset, genesetsUI } = this.props;
    const { genesToAdd } = this.state;

    return (
      <>
        <AnnoDialog
          isActive={genesetsUI.isAddingGenesToGeneset === geneset}
          inputProps={{ "data-testid": `${geneset}:create-label-dialog` }}
          primaryButtonProps={{
            "data-testid": `${geneset}:submit-gene`,
          }}
          title="Add samples to sample set"
          instruction={`Add samples to ${geneset}`}
          cancelTooltipContent="Close this dialog without adding samples to sample set."
          primaryButtonText="Add samples"
          text={genesToAdd}
          validationError={false}
          annoInput={
            <LabelInput
              onChange={this.handleChange}
              inputProps={{
                "data-testid": "add-genes",
                leftIcon: "manually-entered-data",
                intent: "none",
                autoFocus: true,
              }}
              newLabelMessage="New category"
            />
          }
          handleSubmit={this.handleAddGeneToGeneSet}
          handleCancel={this.disableAddGeneMode}
        />
      </>
    );
  }
}

export default AddGeneToGenesetDialogue;
