const Autosave = (
  state = {
    // cell labels
    obsAnnotationSaveInProgress: false,
    lastSavedAnnoMatrix: null,

    // gene sets
    genesetSaveInProgress: false,
    lastSavedGenesets: null,

    // error state
    error: false,
  },
  action: any,
  nextSharedState: any
) => {
  switch (action.type) {
    case "annoMatrix: init complete": {
      return {
        ...state,
        error: false,
        obsAnnotationSaveInProgress: false,
        lastSavedAnnoMatrix: action.annoMatrix,
      };
    }

    case "writable obs annotations - save started": {
      return {
        ...state,
        obsAnnotationSaveInProgress: true,
      };
    }

    case "writable obs annotations - save error": {
      return {
        ...state,
        error: action.message,
        obsAnnotationSaveInProgress: false,
      };
    }

    case "writable obs annotations - save complete": {
      const { lastSavedAnnoMatrix } = action;
      return {
        ...state,
        obsAnnotationSaveInProgress: false,
        error: false,
        lastSavedAnnoMatrix,
      };
    }

    case "geneset: initial load": {
      return {
        ...state,
        genesetSaveInProgress: false,
        lastSavedGenesets: nextSharedState.genesets.genesets,
      };
    }

    case "autosave: genesets started": {
      return {
        ...state,
        genesetSaveInProgress: true,
      };
    }

    case "autosave: genesets error": {
      return {
        ...state,
        genesetSaveInProgress: false,
        error: action.message,
      };
    }

    case "autosave: genesets complete": {
      const { lastSavedGenesets } = action;
      return {
        ...state,
        genesetSaveInProgress: false,
        error: false,
        lastSavedGenesets,
      };
    }

    default:
      return { ...state };
  }
};

export default Autosave;
