const SankeySelection = (
  state={displaySankey: false, categories: {}, sankeyData: null, refresher: false, numChecked: 0},
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
      state.numChecked = numCheckedNow
      return state;
    }
    case "sankey: set": {
      const { category, value } = action;
      state.categories[category] = value;
      state.refresher = !state.refresher;
      const numCheckedNow = Object.values(state.categories).reduce((a, item) => a + item, 0);
      if (numCheckedNow >= 1) {
        state.displaySankey = true;
      } else {
        state.displaySankey = false;
      }
      state.numChecked = numCheckedNow
      return state;
    }    
    case "sankey: set data": {
      const { data } = action;
      state.sankeyData = data;
      return state;
    }
    case "sankey: rename category": {
      const { oldCategoryName, newCategoryName } = action;
      const { categories } = state;
      const newCategories = {};
      delete Object.assign(newCategories, categories, {[newCategoryName]: categories[oldCategoryName] })[oldCategoryName];      
      return {
        ...state,
        categories: newCategories
      };
    }    
    case "sankey: reset": {
      return state={displaySankey: false, categories: {}, sankeyData: null, refresher: false, numChecked: 0};
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
