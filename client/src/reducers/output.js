/*
controller state is not part of the undo/redo history
*/
export const outputController = (
  state = {
    pendingFetch: null,
  },
  action
) => {
  switch (action.type) {
    case "output data: request start": {
      return {
        ...state,
        pendingFetch: action.abortableFetch,
      };
    }
    case "output data: request aborted":
    case "output data: request cancel":
    case "output data: request completed": {
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
export default outputController;