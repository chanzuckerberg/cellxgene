
const SankeySelection = (
  state,
  action
) => {
  switch (action.type) {
    case "sankey: toggle": {
      const { category } = action;
      const value = state?.[category] ?? false;
      state[category] = !value
      console.log(`${category}/${state[category]}`)
      return state;
    }
    default:
      return state;
  }
};

export default SankeySelection;
