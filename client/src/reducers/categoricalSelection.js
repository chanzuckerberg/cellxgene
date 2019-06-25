import { ControlsHelpers } from "../util/stateManager";
import * as globals from "../globals";

function maxCategoryItems(state) {
  return (
    state.config.parameters?.["max-category-items"] ??
    globals.configDefaults.parameters["max-category-items"]
  );
}

const CategoricalSelection = (
  state,
  action,
  nextSharedState,
  prevSharedState
) => {
  switch (action.type) {
    case "initial data load complete (universe exists)":
    case "set World to current selection":
    case "reset World to eq Universe":
    case "set clip quantiles": {
      const { world } = nextSharedState;
      return ControlsHelpers.createCategoricalSelection(
        maxCategoryItems(prevSharedState),
        world
      );
    }

    case "categorical metadata filter select": {
      /*
      Set the specific category in this field to false
      */
      const newCategoryValueSelected = Array.from(
        state[action.metadataField].categoryValueSelected
      );
      newCategoryValueSelected[action.categoryIndex] = true;
      const newCategoricalSelection = {
        ...state,
        [action.metadataField]: {
          ...state[action.metadataField],
          categoryValueSelected: newCategoryValueSelected
        }
      };
      return newCategoricalSelection;
    }

    case "categorical metadata filter deselect": {
      /*
      Set the specific category in this field to false
      */
      const newCategoryValueSelected = Array.from(
        state[action.metadataField].categoryValueSelected
      );
      newCategoryValueSelected[action.categoryIndex] = false;
      const newCategoricalSelection = {
        ...state,
        [action.metadataField]: {
          ...state[action.metadataField],
          categoryValueSelected: newCategoryValueSelected
        }
      };
      return newCategoricalSelection;
    }

    case "categorical metadata filter none of these": {
      /*
      set all categories in this field to false.
      */
      const newCategoricalSelection = {
        ...state,
        [action.metadataField]: {
          ...state[action.metadataField],
          categorySelected: false,
          categoryValueSelected: Array.from(
            state[action.metadataField].categoryValueSelected
          ).fill(false)
        }
      };
      return newCategoricalSelection;
    }

    case "categorical metadata filter all of these": {
      /*
      set all categories in this field to true.
      */
      const newCategoricalSelection = {
        ...state,
        [action.metadataField]: {
          ...state[action.metadataField],
          categorySelected: true,
          categoryValueSelected: Array.from(
            state[action.metadataField].categoryValueSelected
          ).fill(true)
        }
      };
      return newCategoricalSelection;
    }

    default: {
      return state;
    }
  }
};

export default CategoricalSelection;
