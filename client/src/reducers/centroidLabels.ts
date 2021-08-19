const initialState = {
  showLabels: false,
};

const centroidLabels = (
  state = initialState,
  action: any,
  sharedNextState: any
) => {
  const {
    colors: { colorAccessor },
  } = sharedNextState;

  const showLabels = action.showLabels ?? state.showLabels;

  switch (action.type) {
    case "color by categorical metadata":
    case "show centroid labels for category":
      // If colorby is not enabled or labels are not toggled to show
      // then clear the labels and make sure the toggle is off
      return {
        ...state,
        showLabels: colorAccessor && showLabels,
      };

    case "reset centroid labels":
      return initialState;

    default:
      return state;
  }
};

export default centroidLabels;
