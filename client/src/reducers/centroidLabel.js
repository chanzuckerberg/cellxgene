import calcCentroid from "../util/centroid";

/* const initialState = {
  metadataField: "",
  categoryIndex: -1,
  categoryField: "",
  centroidXY: [-1, -1]
}; */

const initialState = {
  labeledCategory: "",
  labels: []
};

const CentroidLabel = (state = initialState, action, sharedNextState) => {
  /* const { categoricalSelection, world, layoutChoice } = sharedNextState;
  const { metadataField, categoryIndex } = action;
  const categoryField =
    categoricalSelection?.[metadataField]?.categoryValues[categoryIndex]; */
  const {
    world,
    layoutChoice,
    categoricalSelection,
    controls
  } = sharedNextState;
  const { metadataField } = action;

  switch (action.type) {
    // case "change graph interaction mode":
    //   if (action.data === "zoom") {
    //     return {
    //       ...initialState,
    //       labeledCategory: "-1"
    //     };
    //   }
    //   return initialState;

    case "show centroid labels for category":
      if (state.labeledCategory === metadataField) {
        return initialState;
      }
      return {
        ...state,
        labeledCategory: metadataField,
        labels: calcCentroid(
          world,
          metadataField,
          layoutChoice.currentDimNames,
          categoricalSelection
        )
      };
    /*     case "category value mouse hover start":
      return {
        ...state,
        metadataField,
        categoryIndex,
        categoryField,
        centroidXY: calcCentroid(  //This function call is computationally heavy and also leading to large GC. Before reimplementation, look into optimization and memoization
          world,
          metadataField,
          categoryField,
          layoutChoice.currentDimNames
        )
      };

    case "category value mouse hover end":
      if (
        metadataField === state.metadataField &&
        categoryIndex === state.categoryIndex
      ) {
        return initialState;
      }
      return state; */

    default:
      return state;
  }
};

export default CentroidLabel;
