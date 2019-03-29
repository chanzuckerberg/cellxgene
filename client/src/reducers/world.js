import _ from "lodash";

import { World, ControlsHelpers } from "../util/stateManager";

const WorldReducer = (
  state = {
    continuousPercentileMin: 0,
    continuousPercentileMax: 100
  },
  action,
  nextSharedState,
  prevSharedState
) => {
  switch (action.type) {
    case "initial data load complete (universe exists)": {
      const { universe } = nextSharedState;
      const world = World.createWorldFromEntireUniverse(
        universe,
        state.continuousPercentileMin,
        state.continuousPercentileMax
      );
      return world;
    }

    case "reset World to eq Universe": {
      return prevSharedState.resetCache.world;
    }

    case "set World to current selection": {
      /* Set viewable world to be the currently selected data */
      const world = World.createWorldFromCurrentSelection(
        action.universe,
        action.world,
        action.crossfilter,
        state.continuousPercentileMin,
        state.continuousPercentileMax
      );
      return world;
    }

    case "expression load success": {
      const { universe } = nextSharedState;
      const universeVarData = universe.varData;
      let worldVarData = state.varData;

      // Load new expression data into the varData dataframes, if
      // not already present.
      _.forEach(action.expressionData, (val, key) => {
        // If not already in world.varData, save sliced expression column
        if (!worldVarData.hasCol(key)) {
          // Slice if world !== universe, else just use whole column.
          // Use the obsAnnotation index as the cut key, as we keep
          // all world dataframes in sync.
          let worldValSlice = val;
          if (!World.worldEqUniverse(state, universe)) {
            worldValSlice = universeVarData
              .subset(state.obsAnnotations.rowIndex.keys(), [key], null)
              .icol(0)
              .asArray();
          }

          // Now build world's varData dataframe
          worldVarData = worldVarData.withCol(
            key,
            worldValSlice,
            state.obsAnnotations.rowIndex
          );
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
      worldVarData = ControlsHelpers.pruneVarDataCache(
        worldVarData,
        allTheGenesWeNeed
      );

      return {
        ...state,
        varData: worldVarData
      };
    }

    default: {
      return state;
    }
  }
};

export default WorldReducer;
