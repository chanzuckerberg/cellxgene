import React from "react";
import { connect } from "react-redux";
import { Button } from "@blueprintjs/core";
import GeneSet from "./geneSet";

import CreateGenesetDialogue from "./menus/createGenesetDialogue";

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

    for (const [name, geneset] of genesets) {
      sets.push(
        <GeneSet
          key={name}
          setGenes={Array.from(geneset.genes.keys())}
          setName={name}
          genesetDescription={geneset.genesetDescription}
        />
      );
    }
    return sets;
  };

  renderDiffexpGeneSets = () => {
    const { differential } = this.props;
    const { diffExp } = differential;
    if (!diffExp) return null;

    // [ [gene, logfoldchange, pval, pval_adj], ...]
    const setGenes = diffExp.map((diffExpGene) => diffExpGene[0]);
    return (
      <GeneSet
        key="Temp DiffExp Set"
        setGenes={setGenes}
        isDiffExp
        diffExp={diffExp}
        setName="Temp DiffExp Set"
      />
    );
  };

  handleActivateCreateGenesetMode = () => {
    const { dispatch } = this.props;
    dispatch({ type: "geneset: activate add new geneset mode" });
  };

  render() {
    return (
      <div>
        <div>
          <div style={{ marginBottom: 10, position: "relative", top: -2 }}>
            <Button
              data-testid="open-create-geneset-dialog"
              onClick={this.handleActivateCreateGenesetMode}
              intent="primary"
            >
              Create new <strong>gene set</strong>
            </Button>
          </div>
          <CreateGenesetDialogue />
        </div>
        <div>{this.renderDiffexpGeneSets()}</div>
        <div>{this.renderGeneSets()}</div>
      </div>
    );
  }
}

export default GeneExpression;
