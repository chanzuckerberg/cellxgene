/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import { Button } from "@blueprintjs/core";
import GeneSet from "./geneSet";

import testGeneSets from "./test_data";
import CreateGenesetDialogue from "./menus/createGenesetDialogue";

@connect((state) => {
  return {
    userDefinedGenes: state.controls.userDefinedGenes,
    differential: state.differential,
    genesets: state.genesets,
  };
})
class GeneExpression extends React.Component {
  renderTestGeneSets = () => {
    const sets = [];

    _.forEach(testGeneSets, (setGenes, setName) => {
      sets.push(
        <GeneSet key={setName} setGenes={setGenes} setName={setName} />
      );
    });

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

    console.log(differential.diffExp);

    return differential.diffExp ? (
      <GeneSet
        key="Temp DiffExp Set"
        setGenes={setGenes}
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
        </div>
        <div>{this.renderDiffexpGeneSets()}</div>
        <div>{this.renderTestGeneSets()}</div>
      </div>
    );
  }
}

export default GeneExpression;
