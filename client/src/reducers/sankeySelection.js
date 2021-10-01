const SankeySelection = (
  state={displaySankey: false, categories: {}, sankeyData: null, refresher: false},
  action
) => {
  switch (action.type) {
    case "sankey: toggle": {
      const { category } = action;
      const value = state?.categories?.[category] ?? false;
      state.categories[category] = !value;
      state.refresher = !state.refresher;
      const numCheckedNow = Object.values(state.categories).reduce((a, item) => a + item, 0);
      if (numCheckedNow >= 1) {
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
    case "sankey: reset": {
      return state={displaySankey: false, categories: {}, sankeyData: null, refresher: false};
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
