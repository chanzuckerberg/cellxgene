const ReembedController = (
	state = {
		pendingFetch: null
	},
	action
) => {
	switch (action.type) {
		case "reembed: request start": {
			return {
				pendingFetch: action.abortableFetch,
				...state
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

export default ReembedController;
