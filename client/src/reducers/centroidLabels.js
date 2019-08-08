import calcCentroid from "../util/centroid";

const initialState = {
  labeledCategory: "",
  labels: []
};

const centroidLabels = (state = initialState, action, sharedNextState) => {
  const { world, layoutChoice, categoricalSelection } = sharedNextState;
  const { metadataField } = action;

  switch (action.type) {
    case "set World to current selection":
      return {
        ...state,
        labels: state.labeledCategory
          ? calcCentroid(
              world,
              state.labeledCategory,
              layoutChoice.currentDimNames,
              categoricalSelection
            )
          : []
      };
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

    default:
      return state;
  }
};

export default centroidLabels;
