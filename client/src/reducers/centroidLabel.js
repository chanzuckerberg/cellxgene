const initialState = {
  metadataField: "",
  categoryIndex: -1
};

const CentroidLabel = (
  state = initialState,
  action,
  nextSharedState,
  prevSharedState
) => {
  const { metadataField, categoryIndex } = action;
  switch (action.type) {
    case "mouse enter":
      console.log(action);
      
      return {
        ...state,
        metadataField,
        categoryIndex
      }

    case "mouse exit":
      console.log(action);
      
        if (metadataField === state.metadataField &&
          categoryIndex === state. categoryIndex) {
          return initialState
        }
        return state;
      default:
      return state
  }
};

export default CentroidLabel;