const initialState = {
  metadataField: "",
  categoryIndex: -1
};

const CentroidLabel = (
  state = initialState,
  action,
  sharedNextState
  ) => {
  const { metadataField, categoryIndex } = action;
  switch (action.type) {
    case "category value mouse hover start":
      console.log(action);
      
      return {
        ...state,
        metadataField,
        categoryIndex
      }

    case "category value mouse hover end":
      console.log(action);
      
        if (metadataField === state.metadataField &&
          categoryIndex === state.categoryIndex) {
          return initialState
        }
        return state;
        
      default:
      return state
  }
};

export default CentroidLabel;