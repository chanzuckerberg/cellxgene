import _ from "lodash";

import { ControlsHelpers } from "../util/stateManager";
import * as globals from "../globals";

function maxCategoryItems(state) {
  return _.get(
    state.config,
    "parameters.max-category-items",
    globals.configDefaults.parameters["max-category-items"]
  );
}

const CategoricalSelectionState = (
  state,
  action,
  nextSharedState,
  prevSharedState
) => {
  switch (action.type) {
    case "initial data load complete (universe exists)":
    case "set World to current selection":
    case "reset World to eq Universe": {
      const { world } = nextSharedState;
      return ControlsHelpers.createCategoricalSelectionState(
        maxCategoryItems(prevSharedState),
        world
      );
    }

    case "categorical metadata filter select": {
      const newCategorySelected = Array.from(
        state[action.metadataField].categorySelected
      );
      newCategorySelected[action.categoryIndex] = true;
      const newCategoricalSelectionState = {
        ...state,
        [action.metadataField]: {
          ...state[action.metadataField],
          categorySelected: newCategorySelected
        }
      };
      return newCategoricalSelectionState;
    }

    case "categorical metadata filter deselect": {
      const newCategorySelected = Array.from(
        state[action.metadataField].categorySelected
      );
      newCategorySelected[action.categoryIndex] = false;
      const newCategoricalSelectionState = {
        ...state,
        [action.metadataField]: {
          ...state[action.metadataField],
          categorySelected: newCategorySelected
        }
      };
      return newCategoricalSelectionState;
    }

    case "categorical metadata filter none of these": {
      const newCategoricalSelectionState = {
        ...state,
        [action.metadataField]: {
          ...state[action.metadataField],
          categorySelected: Array.from(
            state[action.metadataField].categorySelected
          ).fill(false)
        }
      };
      return newCategoricalSelectionState;
    }

    case "categorical metadata filter all of these": {
      const newCategoricalSelectionState = {
        ...state,
        [action.metadataField]: {
          ...state[action.metadataField],
          categorySelected: Array.from(
            state[action.metadataField].categorySelected
          ).fill(true)
        }
      };
      return newCategoricalSelectionState;
    }

    default: {
      return state;
    }
  }
};

export default CategoricalSelectionState;
