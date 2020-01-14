import calcCentroid from "../util/centroid";

const initialState = {
  labels: [],
  toggle: false
};

const centroidLabels = (state = initialState, action, sharedNextState) => {
  const {
    world,
    layoutChoice,
    categoricalSelection,
    colors: { colorAccessor }
  } = sharedNextState;

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
          !!colorAccessor &&
          state.toggle &&
          !!categoricalSelection[colorAccessor]
            ? calcCentroid(
                world.obsAnnotations,
                world.obsLayout,
                colorAccessor,
                layoutChoice.currentDimNames,
                categoricalSelection,
                world.schema.annotations.obsByName
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
          world.obsAnnotations,
          world.obsLayout,
          colorAccessor,
          layoutChoice.currentDimNames,
          categoricalSelection,
          world.schema.annotations.obsByName
        ),
        toggle: action.toggle === undefined ? state.toggle : action.toggle
      };

    case "color by continuous metadata":
      return { ...state, labels: [] };

    case "reset centroid labels":
      return initialState;

    default:
      return state;
  }
};

export default centroidLabels;
