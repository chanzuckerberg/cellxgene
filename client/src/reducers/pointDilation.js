const initialState = {
  metadataField: "",
  categoryField: "",
};

const pointDialation = (state = initialState, action) => {
  const { metadataField, label: categoryField } = action;

  switch (action.type) {
    case "category value mouse hover start":
      return {
        ...state,
        metadataField,
        categoryField,
      };

    case "category value mouse hover end":
      if (
        metadataField === state.metadataField &&
        categoryField === state.categoryField
      ) {
        return initialState;
      }
      return state;

    default:
      return state;
  }
};

export default pointDialation;
