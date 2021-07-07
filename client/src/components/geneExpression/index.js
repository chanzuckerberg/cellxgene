import React from "react";
import { connect } from "react-redux";
import { Button } from "@blueprintjs/core";
import GeneSet from "./geneSet";
import { GenesetHotkeys } from "../hotkeys";

import CreateGenesetDialogue from "./menus/createGenesetDialogue";

@connect((state) => {
  return {
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
          setGenesWithDescriptions={geneset.genes}
          setName={name}
          genesetDescription={geneset.genesetDescription}
        />
      );
    }
    return sets;
  };

  handleActivateCreateGenesetMode = () => {
    const { dispatch } = this.props;
    dispatch({ type: "geneset: activate add new geneset mode" });
  };

  render() {
    const { dispatch, genesets } = this.props;
    return (
      <div>
        <GenesetHotkeys dispatch={dispatch} genesets={genesets} />
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
        <div>{this.renderGeneSets()}</div>
      </div>
    );
  }
}

export default GeneExpression;
