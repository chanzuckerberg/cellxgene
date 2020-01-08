import calcCentroid from "../util/centroid";

const initialState = {
  labels: [],
  toggle: false
};

const centroidLabels = (state = initialState, action, sharedNextState) => {
  const { world, layoutChoice, categoricalSelection, colors } = sharedNextState;
  const { colorAccessor } = colors;

  switch (action.type) {
    case "set World to current selection":
      return {
        ...state,
        labels:
          colorAccessor && !!categoricalSelection[colorAccessor]
            ? calcCentroid(
                world,
                colorAccessor,
                layoutChoice.currentDimNames,
                categoricalSelection
              )
            : []
      };

    case "color by categorical metadata":
    case "show centroid labels for category":
      if (
        !colorAccessor ||
        !(action.toggle || (action.toggle === undefined && state.toggle))
      ) {
        return {
          ...state,
          labels: [],
          toggle: action.toggle === undefined ? state.toggle : action.toggle
        };
      }

      return {
        ...state,
        labels: calcCentroid(
          world,
          colorAccessor,
          layoutChoice.currentDimNames,
          categoricalSelection
        ),
        toggle: action.toggle === undefined ? state.toggle : action.toggle
      };

    case "color by continuous metadata":
      return { ...state, labels: [] };

    default:
      return state;
  }
};

export default centroidLabels;
