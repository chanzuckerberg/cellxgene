/*
controller state is not part of the undo/redo history
*/
export const reembedController = (
  state = {
    pendingFetch: null,
  },
  action
) => {
  switch (action.type) {
    case "reembed: request start": {
      return {
        ...state,
        pendingFetch: action.abortableFetch,
      };
    }
    case "reembed: request aborted":
    case "reembed: request cancel":
    case "reembed: request completed": {
      return {
        ...state,
        pendingFetch: null,
      };
    }
    default: {
      return state;
    }
  }
};

export const reembedParameters = (
  state = {
    dimredParams: {},
    prepParams: {},
    batchParams: {}
  },
  action
) => {
  switch (action.type) {
    case "reembed: set batch correction parameters": {
      const { params } = action;
      return {
        batchParams: params,
      };
    } 
    case "reembed: set dimensionality reduction parameters": {
      const { params } = action;
      return {
        dimredParams: params,
      };
    }     
    case "reembed: set preprocessing parameters": {
      const { params } = action;
      return {
        prepParams: params,
      };
    }     
    default:
      return state;
  }
};