import calcCentroid from "../util/centroid";

const initialState = {
  labels: [],
  showLabels: false,
};

const centroidLabels = (state = initialState, action, sharedNextState) => {
  const {
    world,
    layoutChoice,
    categoricalSelection,
    colors: { colorAccessor },
  } = sharedNextState;

  const showLabels = action.showLabels ?? state.showLabels;

  switch (action.type) {
    case "annotation: label current cell selection":
    case "annotation: label edited":
    case "annotation: delete label":
    case "set layout choice":
    case "set World to current selection":
    case "reset World to eq Universe":
      return {
        ...state,
        labels:
          !!colorAccessor && showLabels && !!categoricalSelection[colorAccessor]
            ? calcCentroid(
                world.obsAnnotations,
                world.obsLayout,
                colorAccessor,
                layoutChoice.currentDimNames,
                categoricalSelection,
                world.schema.annotations.obsByName
              )
            : [],
      };

    case "color by categorical metadata":
    case "show centroid labels for category":
      // If colorby is not enabled or labels are not toggled to show
      // then clear the labels and make sure the toggle is off
      if (!colorAccessor || !showLabels) {
        return {
          ...state,
          labels: [],
          showLabels,
        };
      }

      return {
        ...state,
        labels: calcCentroid(
          world.obsAnnotations,
          world.obsLayout,
          colorAccessor,
          layoutChoice.currentDimNames,
          categoricalSelection,
          world.schema.annotations.obsByName
        ),
        showLabels,
      };

    case "color by continuous metadata":
    case "color by expression":
      return { ...state, labels: [] };

    case "reset centroid labels":
      return initialState;

    default:
      return state;
  }
};

export default centroidLabels;
