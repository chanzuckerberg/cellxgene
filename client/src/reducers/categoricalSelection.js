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
        world,
        ControlsHelpers.selectableCategoryNames(
          world,
          maxCategoryItems(prevSharedState)
        )
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

    case "new user annotation category created": {
      const { world } = nextSharedState;
      const name = action.data;
      return {
        ...state,
        ...ControlsHelpers.createCategoricalSelection(world, [name])
      };
    }

    case "duplicate annotation category": {
      // this code is probably right, but awaiting cleanup of the actions
      // const { world } = nextSharedState;
      // const name = action.metadataField;
      // return {
      //   ...state,
      //   ...ControlsHelpers.createCategoricalSelection(world, [name])
      // };
      return state;
    }

    case "category edited": {
      const name = action.metadataField;
      const newName = action.editedCategoryText;
      const { [name]: catSeln, ...newState } = state;
      newState[newName] = catSeln;
      return newState;
    }

    case "delete category": {
      const name = action.metadataField;
      const { [name]: _, ...newState } = state;
      return newState;
    }

    default: {
      return state;
    }
  }
};

export default CategoricalSelection;
