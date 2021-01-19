/*
Reducers for geneset UI-state.
*/
const GeneSetsUI = (
  state = {
    createGenesetModeActive: false,
    isEditingGenesetName: false,
    genesetAddingNewGenes: null,
  },
  action
) => {
  switch (action.type) {
    case "geneset: activate add new geneset mode": {
      return {
        ...state,
        createGenesetModeActive: true,
      };
    }
    case "geneset: disable create geneset mode": {
      return {
        ...state,
        createGenesetModeActive: false,
      };
    }
    default:
      return state;
  }
};

export default GeneSetsUI;
