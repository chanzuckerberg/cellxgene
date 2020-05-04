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

/*
actual reembedding data is part of the undo/redo history
*/
export const reembedding = (
  state = {
    reembeddings: new Map(),
  },
  action
) => {
  switch (action.type) {
    case "reembed: add reembedding": {
      const { schema, embedding } = action;
      const { name } = schema.name;
      const { reembeddings } = state;
      return {
        ...state,
        reembeddings: new Map(reembeddings).set(name, {
          name,
          schema,
          embedding,
        }),
      };
    }
    case "reembed: clear all reembeddings": {
      return {
        ...state,
        reembeddings: new Map(),
      };
    }
    default: {
      return state;
    }
  }
};
