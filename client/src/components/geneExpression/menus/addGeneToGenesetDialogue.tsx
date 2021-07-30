import React from "react";
import { connect } from "react-redux";
import AnnoDialog from "../../annoDialog";
import LabelInput from "../../labelInput";
import parseBulkGeneString from "../../../util/parseBulkGeneString";
import actions from "../../../actions";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  genesetsUI: (state as any).genesetsUI,
}))
// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
class AddGeneToGenesetDialogue extends React.PureComponent<{}, State> {
  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);
    this.state = {
      genesToAdd: "",
    };
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  disableAddGeneMode = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    dispatch({
      type: "geneset: disable add new genes mode",
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleAddGeneToGeneSet = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'geneset' does not exist on type 'Readonl... Remove this comment to see the full error message
    const { geneset, dispatch } = this.props;
    const { genesToAdd } = this.state;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const genesTmpHardcodedFormat: any = [];
    const genesArrayFromString = parseBulkGeneString(genesToAdd);
    genesArrayFromString.forEach((_gene) => {
      genesTmpHardcodedFormat.push({
        geneSymbol: _gene,
      });
    });
    dispatch(actions.genesetAddGenes(geneset, genesTmpHardcodedFormat));
    dispatch({
      type: "geneset: disable add new genes mode",
    });
    if (e) e.preventDefault();
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleChange = (e: any) => {
    this.setState({ genesToAdd: e });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'geneset' does not exist on type 'Readonl... Remove this comment to see the full error message
    const { geneset, genesetsUI } = this.props;
    const { genesToAdd } = this.state;
    return (
      <>
        <AnnoDialog
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ isActive: boolean; inputProps: { "data-tes... Remove this comment to see the full error message
          isActive={genesetsUI.isAddingGenesToGeneset === geneset}
          inputProps={{ "data-testid": `${geneset}:create-label-dialog` }}
          primaryButtonProps={{
            "data-testid": `${geneset}:submit-gene`,
          }}
          title="Add genes to gene set"
          instruction={`Add genes to ${geneset}`}
          cancelTooltipContent="Close this dialog without adding genes to gene set."
          primaryButtonText="Add genes"
          text={genesToAdd}
          validationError={false}
          annoInput={
            <LabelInput
              // @ts-expect-error ts-migrate(2322) FIXME: Type '{ onChange: (e: any) => void; inputProps: { ... Remove this comment to see the full error message
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
