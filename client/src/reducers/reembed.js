/*
controller state is not part of the undo/redo history
*/
export const reembedController = (
	state = {
		pendingFetch: null
	},
	action
) => {
	switch (action.type) {
		case "reembed: request start": {
			return {
				...state,
				pendingFetch: action.abortableFetch
			};
		}
		case "reembed: request aborted":
		case "reembed: request cancel":
		case "reembed: request completed": {
			return {
				...state,
				pendingFetch: null
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
		reembedding: null
	},
	action
) => {
	switch (action.type) {
		case "reembed: request completed": {
			return {
				...state,
				reembedding: action.reembedding
			};
		}
		default: {
			return state;
		}
	}
};
