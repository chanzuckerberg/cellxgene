import { makeContinuousDimensionName } from "../util/nameCreators";

const ContinuousSelection = (state = {}, action) => {
  switch (action.type) {
    case "reset World to eq Universe": {
      return {};
    }
    case "continuous metadata histogram start":
    case "continuous metadata histogram brush":
    case "continuous metadata histogram end": {
      const name = makeContinuousDimensionName(
        action.continuousNamespace,
        action.selection
      );
      const continuousSelection = { ...state };
      continuousSelection[name] = action.range;
      return continuousSelection;
    }
    default: {
      return state;
    }
  }
};

export default ContinuousSelection;
