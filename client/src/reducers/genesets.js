/*
Reducers for geneset UI-state.
*/
const GeneSets = (
  state = {
    addGenesetModeActive: false,
    isEditingGenesetName: false,
    genesetAddingNewGenes: null,
  },
  action
) => {
  switch (action.type) {
    case "geneset: activate add new geneset mode": {
      return {
        ...state,
        addGenesetModeActive: true,
      };
    }

    default:
      return state;
  }
};

export default GeneSets;
