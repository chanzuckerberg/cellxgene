/*
Reducers for geneset UI-state.
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
    case "geneset: activate add new genes mode": {
      return {
        ...state,
        isAddingGenesToGeneset: action.geneset,
      };
    }
    case "geneset: disable add new genes mode": {
      return {
        ...state,
        isAddingGenesToGeneset: false,
      };
    }
    case "geneset: activate rename geneset mode": {
      return {
        ...state,
        isEditingGenesetName: action.data,
      };
    }
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
