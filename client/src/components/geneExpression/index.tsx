import React from "react";
import { connect } from "react-redux";
import { Button, H4, Icon } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";

import GeneSet from "./geneSet";
import QuickGene from "./quickGene";
import CreateGenesetDialogue from "./menus/createGenesetDialogue";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  genesets: (state as any).genesets.genesets,
}))
// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
class GeneExpression extends React.Component<{}, State> {
  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);
    this.state = { geneSetsExpanded: true };
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  renderGeneSets = () => {
    const sets = [];
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'genesets' does not exist on type 'Readon... Remove this comment to see the full error message
    const { genesets } = this.props;

    for (const [name, geneset] of genesets) {
      sets.push(
        <GeneSet
          key={name}
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ key: any; setGenes: any; setName: any; gen... Remove this comment to see the full error message
          setGenes={geneset.genes}
          setName={name}
          genesetDescription={geneset.genesetDescription}
        />
      );
    }
    return sets;
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleActivateCreateGenesetMode = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    const { geneSetsExpanded } = this.state;
    dispatch({ type: "geneset: activate add new geneset mode" });
    if (!geneSetsExpanded) {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      this.setState((state: any) => ({ ...state, geneSetsExpanded: true }));
    }
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleExpandGeneSets = () => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    this.setState((state: any) => ({
      ...state,
      geneSetsExpanded: !state.geneSetsExpanded,
    }));
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
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
              // @ts-expect-error ts-migrate(2322) FIXME: Type 'string' is not assignable to type 'number | ... Remove this comment to see the full error message
              tabIndex="0"
              data-testclass="geneset-heading-expand"
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
