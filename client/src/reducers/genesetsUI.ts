/*
Reducers for geneset UI-state.
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
const GeneSetsUI = (
  state = {
    createGenesetModeActive: false,
    isEditingGenesetName: false,
    isAddingGenesToGeneset: false,
  },
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  action: any
) => {
  switch (action.type) {
    /**
     * Activate interface for adding a new geneset
     * No params, if the action is fired we flip
     * a boolean here.
     */
    case "geneset: activate add new geneset mode": {
      return {
        ...state,
        createGenesetModeActive: true,
      };
    }
    /**
     * Disable interface for adding a new geneset
     * No params, if the action is fired we flip
     * a boolean here.
     */
    case "geneset: disable create geneset mode": {
      return {
        ...state,
        createGenesetModeActive: false,
      };
    }
    /**
     * Activate the interface for adding new genes to a geneset
     * isAddingGenesToGeneset {
     *  geneset: string, name of geneset
     * },
     */
    case "geneset: activate add new genes mode": {
      return {
        ...state,
        isAddingGenesToGeneset: action.geneset,
      };
    }
    /**
     * Disable the interface for adding new genes to a geneset
     * No params, if the action is fired we flip
     * a boolean here.
     */
    case "geneset: disable add new genes mode": {
      return {
        ...state,
        isAddingGenesToGeneset: false,
      };
    }
    /**
     * Activate the interface for renaming a geneset
     * isEditingGenesetName: {
     *   type: "geneset: activate rename geneset mode",
     *   data: geneset, // a string, name of geneset
     * }
     */
    case "geneset: activate rename geneset mode": {
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
    case "geneset: disable rename geneset mode": {
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
