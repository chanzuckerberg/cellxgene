import { Action } from "redux";
import type { RootState } from ".";
import {
  ContinuousNamespace,
  makeContinuousDimensionName,
} from "../util/nameCreators";

export interface ContinuousSelectionAction extends Action<string> {
  continuousNamespace: ContinuousNamespace;
  selection: string;
  range: [number, number];
}

export interface ContinuousSelectionState {
  [name: string]: [number, number];
}

const ContinuousSelection = (
  state: ContinuousSelectionState = {},
  action: ContinuousSelectionAction
): RootState => {
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
      const { [name]: deletedField, ...newState } = state;
      return newState;
    }
    default: {
      return state;
    }
  }
};

export default ContinuousSelection;
