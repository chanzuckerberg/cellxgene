import _ from "lodash";

import { ControlsHelpers } from "../util/stateManager";

const Universe = (state = null, action, nextSharedState, prevSharedState) => {
  switch (action.type) {
    case "initial data load complete (universe exists)": {
      const { universe } = action;
      return universe;
    }

    case "expression load success": {
      let { varData } = state;

      // Load new expression data into the varData dataframes, if
      // not already present.
      _.forEach(action.expressionData, (val, key) => {
        // If not already in universe.varData, save entire expression column
        if (!varData.hasCol(key)) {
          varData = varData.withCol(key, val);
        }
      });

      // Prune size of varData "cache" if getting out of hand....
      const { userDefinedGenes, diffexpGenes } = prevSharedState;
      const allTheGenesWeNeed = _.uniq(
        [].concat(
          userDefinedGenes,
          diffexpGenes,
          Object.keys(action.expressionData)
        )
      );
      varData = ControlsHelpers.pruneVarDataCache(varData, allTheGenesWeNeed);

      return {
        ...state,
        varData
      };
    }

    default: {
      return state;
    }
  }
};

export default Universe;
