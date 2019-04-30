import { ControlsHelpers } from "../util/stateManager";

const Universe = (state = null, action, nextSharedState, prevSharedState) => {
  switch (action.type) {
    case "initial data load complete (universe exists)": {
      const { universe } = action;
      return universe;
    }

    case "expression load success": {
      let { varData } = state;

      // Lazy load new expression data into the varData dataframe, if
      // not already present.
      //
      Object.entries(action.expressionData).forEach(([key, val]) => {
        // If not already in universe.varData, save entire expression column
        if (!varData.hasCol(key)) {
          varData = varData.withCol(key, val);
        }
      });

      // Prune size of varData "cache" if getting out of hand....
      //
      const { userDefinedGenes, diffexpGenes } = prevSharedState;
      const allTheGenesWeNeed = [
        ...new Set(
          userDefinedGenes,
          diffexpGenes,
          Object.keys(action.expressionData)
        )
      ];
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
