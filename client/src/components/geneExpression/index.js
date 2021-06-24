import React from "react";
import { connect } from "react-redux";
import { Button, H4, Icon } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";

import GeneSet from "./geneSet";
import QuickGene from "./quickGene";
import CreateGenesetDialogue from "./menus/createGenesetDialogue";

@connect((state) => {
  return {
    genesets: state.genesets.genesets,
  };
})
class GeneExpression extends React.Component {
  constructor(props) {
    super(props);
    this.state = { geneSetsExpanded: true };
  }

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
    const { geneSetsExpanded } = this.state;
    dispatch({ type: "geneset: activate add new geneset mode" });
    if (!geneSetsExpanded) {
      this.setState((state) => {
        return { ...state, geneSetsExpanded: true };
      });
    }
  };

  handleExpandGeneSets = () => {
    this.setState((state) => {
      return { ...state, geneSetsExpanded: !state.geneSetsExpanded };
    });
  };

  render() {
    const { geneSetsExpanded } = this.state;
    return (
      <div>
        <QuickGene />
        <div>
          <div
            style={{
              display: "flex",
              flexDirection: "row",
              justifyContent: "space-between",
            }}
          >
            <H4
              role="menuitem"
              tabIndex="0"
              data-testclass="category-expand"
              onKeyPress={this.handleExpandGeneSets}
              style={{
                cursor: "pointer",
              }}
              onClick={this.handleExpandGeneSets}
            >
              Gene Sets{" "}
              {geneSetsExpanded ? (
                <Icon icon={IconNames.CHEVRON_DOWN} />
              ) : (
                <Icon icon={IconNames.CHEVRON_RIGHT} />
              )}
            </H4>

            <div style={{ marginBottom: 10, position: "relative", top: -2 }}>
              <Button
                data-testid="open-create-geneset-dialog"
                onClick={this.handleActivateCreateGenesetMode}
                intent="primary"
              >
                Create new
              </Button>
            </div>
          </div>
          <CreateGenesetDialogue />
        </div>

        {geneSetsExpanded && <div>{this.renderGeneSets()}</div>}
      </div>
    );
  }
}

export default GeneExpression;
