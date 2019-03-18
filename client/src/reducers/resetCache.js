/*
Reducer which caches derived state to be used in a reset
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
