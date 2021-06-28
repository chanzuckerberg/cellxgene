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

export const numPcs = (
  state = {
    npcs: 50,
  },
  action
) => {
  switch (action.type) {
    case "reembed: number of pcs update": {
      const { num } = action;
      return {
        npcs: num,
      };
    }
    case "reembed: reset number of pcs to default":
      return {
        npcs: 50,
      };
    default:
      return state;
  }
};
