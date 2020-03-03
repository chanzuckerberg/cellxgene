import { API } from "../globals";
import { Universe, MatrixFBS } from "../util/stateManager";
import {
	postUserErrorToast,
	postNetworkErrorToast,
	postAsyncSuccessToast,
	postAsyncFailureToast
} from "../components/framework/toasters";

function sleep(ms) {
	return new Promise(resolve => setTimeout(resolve, ms));
}

function abortableFetch(request, opts, timeout = 0) {
	const controller = new AbortController();
	const signal = controller.signal;

	return {
		abort: () => controller.abort(),
		isAborted: () => signal.aborted,
		ready: () => {
			if (timeout) {
				setTimeout(() => controller.abort(), timeout);
			}
			return fetch(request, { ...opts, signal });
		}
	};
}

async function doReembedFetch(dispatch, getState) {
	const state = getState();
	let cells = state.world.obsAnnotations.rowIndex.keys();

	// These lines ensure that we convert any TypedArray to an Array.
	// This is necessary because JSON.stringify() does some very strange
	// things with TypedArrays (they are marshalled to JSON objects, rather
	// than being marshalled as a JSON array).
	cells = Array.isArray(cells) ? cells : Array.from(cells);

	const af = abortableFetch(
		`${API.prefix}${API.version}layout/obs`,
		{
			method: "PUT",
			headers: new Headers({
				Accept: "application/octet-stream",
				"Content-Type": "application/json"
			}),
			body: JSON.stringify({
				method: "umap",
				filter: { obs: { index: cells } }
			}),
			credentials: "include"
		},
		60000 // 1 minute timeout
	);
	dispatch({
		type: "reembed: request start",
		abortableFetch: af
	});
	const res = await af.ready();

	if (
		res.ok &&
		res.headers.get("Content-Type").includes("application/octet-stream")
	) {
		return res;
	}

	// else an error
	let msg = `Unexpected HTTP response ${res.status}, ${res.statusText}`;
	const body = await res.text();
	if (body && body.length > 0) {
		msg = `${msg} -- ${body}`;
	}
	postNetworkErrorToast(msg);
	throw new Error(msg);
}

/*
below are dispatch-able
*/

function cancelReembed() {
	return (dispatch, getState) => {
		const state = getState();
		const { pendingFetch } = state.reembedController;
		dispatch({ type: "reembed: request cancel" });
		pendingFetch?.abort();
	};
}

export function requestReembed() {
	return async (dispatch, getState) => {
		try {
			const res = await doReembedFetch(dispatch, getState);
			const buffer = await res.arrayBuffer();
			const df = Universe.matrixFBSToDataframe(buffer);
			dispatch({
				type: "reembed: request completed",
				dataframe: df
			});
			postAsyncSuccessToast("Re-embedding has completed.");
		} catch (error) {
			dispatch({
				type: "reembed: request aborted"
			});
			if (error.name === "AbortError") {
				postAsyncFailureToast("Re-embedding: took too long and was aborted.");
			} else {
				postNetworkErrorToast(`Re-embedding: ${error.message}`);
			}
			console.log("Reembed exception:", error, error.name, error.message);
		}
	};
}

requestReembed.cancel = cancelReembed;
