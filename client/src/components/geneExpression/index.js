import React from "react";
import { connect } from "react-redux";
import { Button, Classes } from "@blueprintjs/core";
import GeneSet from "./geneSet";

import CreateGenesetDialogue from "./menus/createGenesetDialogue";

@connect((state) => {
  return {
    genesets: state.genesets.genesets,
    skeleton: state.skeleton.skeleton,
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
    const { skeleton } = this.props;
    return (
      <div>
        <div>
          <div style={{ marginBottom: 10, position: "relative", top: -2 }}>
            <Button
              data-testid="open-create-geneset-dialog"
              onClick={this.handleActivateCreateGenesetMode}
              intent="primary"
              className={skeleton ? Classes.SKELETON : null}
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
