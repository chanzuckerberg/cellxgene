const SankeySelection = (
  state={displaySankey: false, categories: {}, sankeyData: null},
  action
) => {
  switch (action.type) {
    case "sankey: toggle": {
      const { category } = action;
      const numChecked = Object.values(state.categories).reduce((a, item) => a + item, 0);
      const value = state?.categories?.[category] ?? false;
      if (numChecked > 1 && !value){
        return state;
      }
      state.categories[category] = !value;
      const numCheckedNow = Object.values(state.categories).reduce((a, item) => a + item, 0);
      if (numCheckedNow == 2) {
        state.displaySankey = true;
      } else {
        state.displaySankey = false;
      }
      return state;
    }
    case "sankey: set data": {
      const { data } = action;
      state.sankeyData = data;
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
