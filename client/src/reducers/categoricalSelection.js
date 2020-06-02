import { ControlsHelpers as CH } from "../util/stateManager";

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
      const newState = CH.createCategoricalSelection(
        world,
        CH.selectableCategoryNames(
          world.schema,
          CH.maxCategoryItems(prevSharedState.config)
        )
      );
      return newState;
    }

    case "universe: column load success": {
      const { dim } = action;
      if (dim !== "obsAnnotations") return state;

      const { dataframe } = action;
      const { world } = nextSharedState;
      const names = CH.selectableCategoryNames(
        world.schema,
        CH.maxCategoryItems(prevSharedState.config),
        dataframe.colIndex.labels()
      );
      if (names.length === 0) return state;
      return {
        ...state,
        ...CH.createCategoricalSelection(world, names),
      };
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
          categoryValueSelected: newCategoryValueSelected,
        },
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
          categoryValueSelected: newCategoryValueSelected,
        },
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
          ).fill(false),
        },
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
          ).fill(true),
        },
      };
      return newCategoricalSelection;
    }

    case "annotation: create category": {
      const { world } = nextSharedState;
      const name = action.data;
      return {
        ...state,
        ...CH.createCategoricalSelection(world, [name]),
      };
    }

    case "annotation: category edited": {
      const name = action.metadataField;
      const newName = action.newCategoryText;
      const { [name]: catSeln, ...newState } = state;
      newState[newName] = catSeln;
      return newState;
    }

    case "annotation: delete category": {
      const name = action.metadataField;
      const { [name]: _, ...newState } = state;
      return newState;
    }

    case "annotation: label current cell selection":
    case "annotation: add new label to category":
    case "annotation: label edited":
    case "annotation: delete label": {
      /* need to rebuild the state for this annotation */
      const { world } = nextSharedState;
      const name = action.metadataField;
      const { [name]: _, ...partialState } = state;
      return {
        ...partialState,
        ...CH.createCategoricalSelection(world, [name]),
      };
    }

    default: {
      return state;
    }
  }
};

export default CategoricalSelection;
