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
      genesToAdd: "",
    };
  }

  disableAddGeneMode = () => {
    const { dispatch } = this.props;
    dispatch({
      type: "geneset: disable add new genes mode",
    });
  };

  handleAddGeneToGeneSet = (e) => {
    const { geneset, dispatch } = this.props;
    const { genesToAdd } = this.state;
    dispatch({
      type: "geneset: add genes",
      name: geneset,
      genes: genesToAdd.split(","),
    });
    dispatch({
      type: "geneset: disable add new genes mode",
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
            "data-testid": `${geneset}:submit-label`,
          }}
          title="Add genes to geneset"
          instruction={`Add gene to ${geneset}`}
          cancelTooltipContent="Close this dialog without adding genes to geneset."
          primaryButtonText="Add genes"
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
