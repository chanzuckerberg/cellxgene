const Autosave = (
	state = {
		saveInProgress: false,
		error: false,
		lastSavedObsAnnotations: null
	},
	action,
	nextSharedState
) => {
	switch (action.type) {
		case "initial data load complete (universe exists)": {
			/* don't save on init */
			const { universe } = nextSharedState;
			return {
				...state,
				error: false,
				saveInProgress: false,
				lastSavedObsAnnotations: universe.obsAnnotations
			};
		}

		case "writable obs annotations - save started": {
			return {
				...state,
				saveInProgress: true
			};
		}

		case "writable obs annotations - save error": {
			const { message } = action;
			return {
				...state,
				error: message,
				saveInProgress: false
			};
		}

		case "writable obs annotations - save complete": {
			const lastSavedObsAnnotations = action.obsAnnotations;
			return {
				...state,
				saveInProgress: false,
				error: false,
				lastSavedObsAnnotations
			};
		}

		default:
			return { ...state };
	}
};

export default Autosave;
