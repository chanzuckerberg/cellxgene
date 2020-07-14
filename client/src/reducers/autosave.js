const Autosave = (
  state = {
    saveInProgress: false,
    error: false,
    lastSavedAnnoMatrix: null,
  },
  action
) => {
  switch (action.type) {
    case "annoMatrix: init complete": {
      return {
        ...state,
        error: false,
        saveInProgress: false,
        lastSavedAnnoMatrix: action.annoMatrix,
      };
    }

    case "writable obs annotations - save started": {
      return {
        ...state,
        saveInProgress: true,
      };
    }

    case "writable obs annotations - save error": {
      return {
        ...state,
        error: action.message,
        saveInProgress: false,
      };
    }

    case "writable obs annotations - save complete": {
      const { lastSavedAnnoMatrix } = action;
      return {
        ...state,
        saveInProgress: false,
        error: false,
        lastSavedAnnoMatrix,
      };
    }

    default:
      return { ...state };
  }
};

export default Autosave;
