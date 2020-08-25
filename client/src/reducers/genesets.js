/*
Reducers for geneset UI-state.
*/
const GeneSets = (
  state = {
    createGenesetModeActive: false,
    isEditingGenesetName: false,
    genesetAddingNewGenes: null,
    userCreatedGenesets: {},
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
    case "geneset: create": {
      return {
        ...state,
        userCreatedGenesets: Object.assign(
          {},
          { [action.genesetName]: {} },
          ...state.userCreatedGenesets
        ),
      };
    }
    default:
      return state;
  }
};

export default GeneSets;
