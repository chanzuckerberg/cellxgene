const Metadata = (
  state = {
    singleValueCategories: new Map(),
  },
  action
) => {
  const { category, value } = action;
  const { singleValueCategories } = state;
  switch (action.type) {
    case "add singleValueCategory":
      return {
        ...state,
        singleValueCategories: singleValueCategories.set(category, value),
      };
    default:
      return state;
  }
};

export default Metadata;
