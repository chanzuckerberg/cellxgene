/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import { connect } from "react-redux";
import { Button } from "@blueprintjs/core";
import GeneSet from "./geneSet";

import CreateGenesetDialogue from "./menus/createGenesetDialogue";
import EditGenesetNameDialogue from "./menus/editGenesetNameDialogue";

@connect((state) => {
  return {
    differential: state.differential,
    genesets: state.genesets.genesets,
  };
})
class GeneExpression extends React.Component {
  renderGeneSets = () => {
    const sets = [];
    const { genesets } = this.props;

    for (const [name, genes] of genesets) {
      sets.push(<GeneSet key={name} setGenes={[...genes]} setName={name} />);
    }
    return sets;
  };

  renderDiffexpGeneSets = () => {
    const { differential } = this.props;

    const setGenes = [];

    if (differential.diffExp) {
      differential.diffExp.forEach((diffexpGene) => {
        setGenes.push(diffexpGene[0]);
      });
    }

    return differential.diffExp ? (
      <GeneSet
        key="Temp DiffExp Set"
        setGenes={setGenes}
        isDiffexp
        setName="Temp DiffExp Set"
      />
    ) : null;
  };

  handleActivateCreateGenesetMode = () => {
    const { dispatch } = this.props;
    dispatch({ type: "geneset: activate add new geneset mode" });
  };

  render() {
    const geneSetsFeatureEnabledTODO = true;
    return (
      <div>
        <div>
          {geneSetsFeatureEnabledTODO ? (
            <div style={{ marginBottom: 10, position: "relative", top: -2 }}>
              <Button
                data-testid="open-create-geneset-dialog"
                onClick={this.handleActivateCreateGenesetMode}
                intent="primary"
              >
                Create new <strong>gene set</strong>
              </Button>
            </div>
          ) : null}
          <CreateGenesetDialogue />
          <EditGenesetNameDialogue />
        </div>
        <div>{this.renderDiffexpGeneSets()}</div>
        <div>{this.renderGeneSets()}</div>
      </div>
    );
  }
}

export default GeneExpression;
