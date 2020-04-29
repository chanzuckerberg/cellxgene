/*
Reducer which caches derived state to be used in a reset or other 
recomputation.  Add stuff here you want stashed at init time (or whenever),
for later use.

Currently this only caches the baseline (full universe) crossfilter, 
which improves Reset UI performance.
*/
const ResetCacheReducer = (
  state = {
    crossfilter: null,
  },
  action,
  nextSharedState
) => {
  switch (action.type) {
    case "initial data load complete (universe exists)": {
      const { crossfilter } = nextSharedState;
      return {
        ...state,
        crossfilter,
      };
    }
    default: {
      return state;
    }
  }
};

export default ResetCacheReducer;
