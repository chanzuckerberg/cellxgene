const initialState = {
  metadataField: "",
  categoryField: "",
};

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
const pointDialation = (state = initialState, action: any) => {
  const { metadataField, label: categoryField } = action;

  switch (action.type) {
    case "category value mouse hover start":
      return {
        ...state,
        metadataField,
        categoryField,
      };

    case "category value mouse hover end":
      if (
        metadataField === state.metadataField &&
        categoryField === state.categoryField
      ) {
        return initialState;
      }
      return state;

    default:
      return state;
  }
};

export default pointDialation;
