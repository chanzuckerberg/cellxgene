const initialState = {
  metadataField: "",
  categoryField: ""
};

const pointDialation = (state = initialState, action, sharedNextState) => {
  const { categoricalSelection } = sharedNextState;
  const { metadataField, categoryIndex } = action;
  const categoryField =
    action.categoryField ||
    categoricalSelection?.[metadataField]?.categoryValues[categoryIndex];

  switch (action.type) {
    case "category value mouse hover start":
      return {
        ...state,
        metadataField,
        categoryField
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
