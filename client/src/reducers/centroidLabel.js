import calcCentroid from "../util/centroid";

const initialState = {
  metadataField: 0,
  categoryField: 0,
  centroidXY: [-1, -1]
};

const CentroidLabel = (state = initialState, action, sharedNextState) => {
  const { categoricalSelection, world, layoutChoice } = sharedNextState;
  const { metadataField, categoryIndex } = action;
  // TODO: Use the es7 syntatic sugar
  const categoryField =
    categoricalSelection && categoricalSelection[metadataField]
      ? categoricalSelection[metadataField].categoryValues[categoryIndex]
      : "";
  switch (action.type) {
    case "category value mouse hover start":
      console.log(action);

      return {
        ...state,
        metadataField,
        categoryIndex,
        centroidXY: calcCentroid(
          world,
          metadataField,
          categoryField,
          layoutChoice.currentDimNames
        )
      };

    case "category value mouse hover end":
      console.log(action);

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
