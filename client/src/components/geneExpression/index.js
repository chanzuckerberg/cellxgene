// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import HistogramBrush from "../brushableHistogram";
import * as globals from "../../globals";
import AddGenes from "./addGenes";

@connect((state) => {
  return {
    obsAnnotations: state.world?.obsAnnotations,
    userDefinedGenes: state.controls.userDefinedGenes,
    userDefinedGenesLoading: state.controls.userDefinedGenesLoading,
    world: state.world,
    colorAccessor: state.colors.colorAccessor,
    differential: state.differential,
  };
})
class GeneExpression extends React.Component {
  render() {
    const { world, userDefinedGenes, differential } = this.props;
    const varIndexName = world?.schema?.annotations?.var?.index;
    const varIndex = world?.varAnnotations?.col(varIndexName)?.asArray();

    // may still be loading!
    if (!varIndex) return null;

    return (
      <div
        style={{
          borderBottom: `1px solid ${globals.lighterGrey}`,
        }}
      >
        <div>
          <AddGenes />
          {world && userDefinedGenes.length > 0
            ? _.map(userDefinedGenes, (geneName, index) => {
                const values = world.varData.col(geneName);
                if (!values) {
                  return null;
                }
                const summary = values.summarize();
                return (
                  <HistogramBrush
                    key={geneName}
                    field={geneName}
                    zebra={index % 2 === 0}
                    ranges={summary}
                    isUserDefined
                  />
                );
              })
            : null}
        </div>
        <div>
          {differential.diffExp
            ? _.map(differential.diffExp, (value, index) => {
                const name = world.varAnnotations.at(value[0], varIndexName);
                const values = world.varData.col(name);
                if (!values) {
                  return null;
                }
                const summary = values.summarize();
                return (
                  <HistogramBrush
                    key={name}
                    field={name}
                    zebra={index % 2 === 0}
                    ranges={summary}
                    isDiffExp
                    logFoldChange={value[1]}
                    pval={value[2]}
                    pvalAdj={value[3]}
                  />
                );
              })
            : null}
        </div>
      </div>
    );
  }
}

export default GeneExpression;
