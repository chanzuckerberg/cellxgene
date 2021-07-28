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
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
const CategoricalSelection = (
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  state: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  action: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  nextSharedState: any
) => {
  switch (action.type) {
    case "initial data load complete":
    case "subset to selection":
    case "reset subset":
    case "set clip quantiles": {
      const { annoMatrix } = nextSharedState;
      const newState = CH.createCategoricalSelection(
        // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
        CH.selectableCategoryNames(annoMatrix.schema)
      );
      return newState;
    }

    case "categorical metadata filter select":
    case "categorical metadata filter deselect":
    case "categorical metadata filter none of these":
    case "categorical metadata filter all of these": {
      const { metadataField, labelSelectionState } = action;
      return {
        ...state,
        [metadataField]: labelSelectionState,
      };
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
