const Leiden = (
  state={res: 1.0},
  action
) => {
  switch (action.type) {
    case "leiden: set resolution": {
      const { res } = action;
      state.res = res;
      return state;
    }
    default:
      return state;
  }
};

export const leidenController = (
  state = {
    pendingFetch: null,
  },
  action
) => {
  switch (action.type) {
    case "leiden: request start": {
      return {
        ...state,
        pendingFetch: action.abortableFetch,
      };
    }
    case "leiden: request aborted":
    case "leiden: request cancel":
    case "leiden: request completed": {
      return {
        ...state,
        pendingFetch: null,
      };
    }
    default: {
      return state;
    }
  }
};

export default Leiden;
