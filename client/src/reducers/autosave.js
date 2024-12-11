const Autosave = (
  state = {
    // protein labels
    obsAnnotationSaveInProgress: false,
    lastSavedAnnoMatrix: null,

    // sample sets
    genesetSaveInProgress: false,
    lastSavedGenesets: null,

    // error state
    error: false,
  },
  action,
  nextSharedState
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

    case "sampleset: initial load": {
      return {
        ...state,
        genesetSaveInProgress: false,
        lastSavedGenesets: nextSharedState.genesets.genesets,
      };
    }

    case "autosave: samplesets started": {
      return {
        ...state,
        genesetSaveInProgress: true,
      };
    }

    case "autosave: samplesets error": {
      return {
        ...state,
        genesetSaveInProgress: false,
        error: action.message,
      };
    }

    case "autosave: samplesets complete": {
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
