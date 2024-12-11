/*
Reducers for sampleset UI-state.
*/
const GeneSetsUI = (
  state = {
    createGenesetModeActive: false,
    isEditingGenesetName: false,
    isAddingGenesToGeneset: false,
  },
  action
) => {
  switch (action.type) {
    /**
     * Activate interface for adding a new sampleset
     * No params, if the action is fired we flip
     * a boolean here.
     */
    case "sampleset: activate add new sampleset mode": {
      return {
        ...state,
        createGenesetModeActive: true,
      };
    }
    /**
     * Disable interface for adding a new sampleset
     * No params, if the action is fired we flip
     * a boolean here.
     */
    case "sampleset: disable create sampleset mode": {
      return {
        ...state,
        createGenesetModeActive: false,
      };
    }
    /**
     * Activate the interface for adding new samples to a sampleset
     * isAddingGenesToGeneset {
     *  sampleset: string, name of geneset
     * },
     */
    case "sampleset: activate add new samples mode": {
      return {
        ...state,
        isAddingGenesToGeneset: action.geneset,
      };
    }
    /**
     * Disable the interface for adding new samples to a sampleset
     * No params, if the action is fired we flip
     * a boolean here.
     */
    case "sampleset: disable add new samples mode": {
      return {
        ...state,
        isAddingGenesToGeneset: false,
      };
    }
    /**
     * Activate the interface for renaming a sampleset
     * isEditingGenesetName: {
     *   type: "sampleset: activate rename sampleset mode",
     *   data: geneset, // a string, name of sampleset
     * }
     */
    case "sampleset: activate rename sampleset mode": {
      return {
        ...state,
        isEditingGenesetName: action.data,
      };
    }
    /**
     * Disable the interface for renaming a geneset
     * No params, if the action is fired we flip
     * a boolean here.
     */
    case "sampleset: disable rename sampleset mode": {
      return {
        ...state,
        isEditingGenesetName: false,
      };
    }

    default:
      return state;
  }
};

export default GeneSetsUI;
