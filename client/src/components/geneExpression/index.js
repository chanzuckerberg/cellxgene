/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import HistogramBrush from "../brushableHistogram";
import * as globals from "../../globals";
import AddGenes from "./addGenes";

@connect((state) => {
  return {
    userDefinedGenes: state.controls.userDefinedGenes,
    differential: state.differential,
  };
})
class GeneExpression extends React.Component {
  render() {
    const { userDefinedGenes, differential } = this.props;
    return (
      <div
        style={{
          borderBottom: `1px solid ${globals.lighterGrey}`,
        }}
      >
        <div>
          <AddGenes />
          {userDefinedGenes.length > 0
            ? _.map(userDefinedGenes, (geneName, index) => {
                return (
                  <HistogramBrush
                    key={geneName}
                    field={geneName}
                    zebra={index % 2 === 0}
                    isUserDefined
                  />
                );
              })
            : null}
        </div>
        <div>
          {differential.diffExp
            ? _.map(differential.diffExp, (value, index) => {
                return (
                  <HistogramBrush
                    key={value[0]}
                    field={value[0]}
                    zebra={index % 2 === 0}
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
