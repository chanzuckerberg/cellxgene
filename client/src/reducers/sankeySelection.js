
const SankeySelection = (
  state={displaySankey: false, categories: {}},
  action
) => {
  switch (action.type) {
    case "sankey: toggle": {
      const { category } = action;
      const value = state?.categories?.[category] ?? false;
      state.categories[category] = !value;
      const numChecked = Object.values(state.categories).reduce((a, item) => a + item, 0);
      if (numChecked > 1) {
        state.displaySankey = true;
      } else {
        state.displaySankey = false;
      }
      return state;
    }
    default:
      return state;
  }
};

export const sankeyController = (
  state = {
    pendingFetch: null,
  },
  action
) => {
  switch (action.type) {
    case "sankey: request start": {
      return {
        ...state,
        pendingFetch: action.abortableFetch,
      };
    }
    case "sankey: request aborted":
    case "sankey: request cancel":
    case "sankey: request completed": {
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

export default SankeySelection;
