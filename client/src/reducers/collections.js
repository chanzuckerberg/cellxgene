const Collections = (
  /*
  collections reducer, modifies portal collections-related state. 
   */
  state = {
    // data loading flag
    loading: true,
    error: null,

    collection: null,
    selectedDatasetId: null,
  },
  action
) => {
  switch (action.type) {
    case "initial data load start":
      return {
        ...state,
        loading: true,
        error: null,
      };
    case "collection load complete":
      return {
        ...state,
        loading: false,
        error: null,
        collection: action.collection,
        selectedDatasetId: action.selectedDatasetId,
      };
    case "initial data load error":
      return {
        ...state,
        error: action.error,
      };
    default:
      return state;
  }
};

export default Collections;
