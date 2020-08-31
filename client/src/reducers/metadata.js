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
      if (value)
        return {
          ...state,
          singleValueCategories: singleValueCategories.set(category, value),
        };
      return state;
    default:
      return state;
  }
};

export default Metadata;
