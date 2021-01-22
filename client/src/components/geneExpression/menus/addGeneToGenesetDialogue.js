import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../categorical/annoDialog";
import LabelInput from "../../categorical/labelInput";

@connect((state) => ({
  genesetsUI: state.genesetsUI,
}))
class AddGeneToGenesetDialogue extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      //   genesToAdd: [],
    };
  }

  disableAddGeneMode = () => {
    const { dispatch } = this.props;
    dispatch({
      type: "geneset: disable add new genes mode",
    });
  };

  render() {
    const { geneset, genesetsUI } = this.props;

    return (
      <>
        <AnnoDialog
          isActive={genesetsUI.isAddingGenesToGeneset === geneset}
          inputProps={{ "data-testid": `${geneset}:create-label-dialog` }}
          primaryButtonProps={{
            "data-testid": `${geneset}:submit-label`,
          }}
          title="Add gene to geneset"
          instruction={`Add gene to ${geneset}`}
          cancelTooltipContent="Close this dialog without adding genes to geneset."
          primaryButtonText="Add Gene(s)"
          text="todo"
          validationError={() => {
            /* todo genesets this.labelNameError(genesetName) */
            return false;
          }}
          annoInput={
            <LabelInput
              onChange={this.handleChange}
              onSelect={this.handleSelect}
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
