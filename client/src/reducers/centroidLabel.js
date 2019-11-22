const initialState = {
  metadataField: "",
  categoryIndex: -1,
  categoryField: "",
  centroidXY: [-1, -1]
};

const CentroidLabel = (state = initialState, action, sharedNextState) => {
  const { categoricalSelection } = sharedNextState;
  const { metadataField, categoryIndex } = action;
  const categoryField =
    categoricalSelection?.[metadataField]?.categoryValues[categoryIndex];
  switch (action.type) {
    case "category value mouse hover start":
      return {
        ...state,
        metadataField,
        categoryIndex,
        categoryField,
        centroidXY: null
      };

    case "category value mouse hover end":
      if (
        metadataField === state.metadataField &&
        categoryIndex === state.categoryIndex
      ) {
        return initialState;
      }
      return state;

    default:
      return state;
  }
};

export default CentroidLabel;
