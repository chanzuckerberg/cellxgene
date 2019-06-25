/*
Reducer which caches derived state to be used in a reset or other 
recomputation.

Currently this only caches the baseline (full universe) world & crossfilter, 
for use in a Reset.
*/
const ResetCacheReducer = (
  state = {
    world: null,
    crossfilter: null
  },
  action,
  nextSharedState
) => {
  switch (action.type) {
    case "initial data load complete (universe exists)": {
      const { world, crossfilter } = nextSharedState;
      return {
        ...state,
        world,
        crossfilter
      };
    }
    default: {
      return state;
    }
  }
};

export default ResetCacheReducer;
