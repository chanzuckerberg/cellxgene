import { makeContinuousDimensionName } from "../util/nameCreators";

const ContinuousSelection = (state = {}, action: any) => {
  switch (action.type) {
    case "reset subset":
    case "subset to selection":
    case "set clip quantiles": {
      return {};
    }
    case "continuous metadata histogram start":
    case "continuous metadata histogram brush":
    case "continuous metadata histogram end": {
      const name = makeContinuousDimensionName(
        action.continuousNamespace,
        action.selection
      );
      return {
        ...state,
        [name]: action.range,
      };
    }
    case "continuous metadata histogram cancel": {
      const name = makeContinuousDimensionName(
        action.continuousNamespace,
        action.selection
      );
      // @ts-expect-error ts-migrate(2537) FIXME: Type '{}' has no matching index signature for type... Remove this comment to see the full error message
      const { [name]: deletedField, ...newState } = state;
      return newState;
    }
    default: {
      return state;
    }
  }
};

export default ContinuousSelection;
