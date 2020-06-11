import { ControlsHelpers as CH } from "../util/stateManager";

/*
State is an object, with a key for each categorical annotation, and a
value which is a Map of label->t/f, reflecting selection state for the label.

Label state default (if missing) is up to the component, but typically true.

{
  "louvain": Map(),
  ...
}
*/
const CategoricalSelection = (state, action, nextSharedState) => {
  switch (action.type) {
    case "initial data load complete (universe exists)":
    case "set World to current selection":
    case "reset World to eq Universe":
    case "set clip quantiles": {
      const { world } = nextSharedState;
      const newState = CH.createCategoricalSelection(
        CH.selectableCategoryNames(world.schema)
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
        dataframe.colIndex.labels()
      );
      if (names.length === 0) return state;
      return {
        ...state,
        ...CH.createCategoricalSelection(names),
      };
    }

    case "categorical metadata filter select": {
      /*
      Set the specific category in this field to false
      */
      const { metadataField, label } = action;
      const newSelected = new Map(state[metadataField]);
      newSelected.set(label, true);
      const newCategoricalSelection = {
        ...state,
        [action.metadataField]: newSelected,
      };
      return newCategoricalSelection;
    }

    case "categorical metadata filter deselect": {
      /*
      Set the specific category in this field to false
      */
      const { metadataField, label } = action;
      const newSelected = new Map(state[metadataField]);
      newSelected.set(label, false);
      const newCategoricalSelection = {
        ...state,
        [action.metadataField]: newSelected,
      };
      return newCategoricalSelection;
    }

    case "categorical metadata filter none of these": {
      /*
      set all categories in this field to false.
      */
      const { metadataField, labels } = action;
      const { selected } = state[metadataField];
      const newSelected = new Map(selected);
      labels.forEach((label) => newSelected.set(label, false));
      const newCategoricalSelection = {
        ...state,
        [action.metadataField]: newSelected,
      };
      return newCategoricalSelection;
    }

    case "categorical metadata filter all of these": {
      /*
      set all categories in this field to true.
      */
      const { metadataField, labels } = action;
      const { selected } = state[metadataField];
      const newSelected = new Map(selected);
      labels.forEach((label) => newSelected.set(label, true));
      const newCategoricalSelection = {
        ...state,
        [action.metadataField]: newSelected,
      };
      return newCategoricalSelection;
    }

    case "annotation: create category": {
      const name = action.data;
      return {
        ...state,
        ...CH.createCategoricalSelection([name]),
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
      const name = action.metadataField;
      const { [name]: _, ...partialState } = state;
      return {
        ...partialState,
        ...CH.createCategoricalSelection([name]),
      };
    }

    default: {
      return state;
    }
  }
};

export default CategoricalSelection;
