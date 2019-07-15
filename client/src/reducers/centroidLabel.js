import calcCentroid from "../util/centroid";

const initialState = {
  metadataField: "",
  categoryIndex: -1,
  categoryField: "",
  centroidXY: [-1, -1]
};

const CentroidLabel = (state = initialState, action, sharedNextState) => {
  const { categoricalSelection, world, layoutChoice } = sharedNextState;
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
        centroidXY: null /* calcCentroid(  This function call is computationally heavy and also leading to large GC. Before reimplementation, look into optimization and memoization
          world,
          metadataField,
          categoryField,
          layoutChoice.currentDimNames
        ) */
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
