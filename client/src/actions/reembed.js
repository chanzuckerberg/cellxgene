import { API } from "../globals";

import {
  postNetworkErrorToast,
  postAsyncSuccessToast,
  postAsyncFailureToast,
} from "../components/framework/toasters";
import { _switchEmbedding } from "./embedding";
import { subsetAction } from "./viewStack";
import { AnnoMatrixLoader } from "../annoMatrix";


function abortableFetch(request, opts, timeout = 0) {
  const controller = new AbortController();
  const { signal } = controller;

  return {
    abort: () => controller.abort(),
    isAborted: () => signal.aborted,
    ready: () => {
      if (timeout) {
        setTimeout(() => controller.abort(), timeout);
      }
      return fetch(request, { ...opts, signal });
    },
  };
}

async function doReembedFetch(dispatch, getState, reembedParams,parentName,embName) {
  const state = getState();
  let cells = state.annoMatrix.rowIndex.labels();
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
        "Content-Type": "application/json",
      }),
      body: JSON.stringify({
        method: "umap",
        filter: { obs: { index: cells } },
        params: reembedParams,
        parentName: parentName,
        embName: embName,
      }),
      credentials: "include",
    },
    6000000 // 1 minute timeout
  );
  dispatch({
    type: "reembed: request start",
    abortableFetch: af,
  });
  const res = await af.ready();

  if (res.ok && res.headers.get("Content-Type").includes("application/json")) {
    // #TOOD: TRIGGER CELL SUBSETTING AND AWAIT RESULTS!
    return res;
  }

  // else an error
  let msg = `Unexpected HTTP response ${res.status}, ${res.statusText}`;
  const body = await res.text();
  if (body && body.length > 0) {
    msg = `${msg} -- ${body}`;
  }
  throw new Error(msg);
}
/*
functions below are dispatch-able
*/
export function requestReembed(reembedParams,parentName,embName) {
  return async (dispatch, getState) => {
    try {
      await dispatch(subsetAction());
      const res = await doReembedFetch(dispatch, getState, reembedParams,parentName,embName);
      const schema = await res.json();
      dispatch({
        type: "reembed: request completed",
      });

      const {
        annoMatrix: prevAnnoMatrix,
        obsCrossfilter: prevCrossfilter,
        layoutChoice,
      } = getState();
      const base = prevAnnoMatrix.base().addEmbedding(schema);

      const [annoMatrix, obsCrossfilter] = await _switchEmbedding(
        base,
        prevCrossfilter,
        layoutChoice.current,
        schema.name
      );
      dispatch({
        type: "reembed: add reembedding",
        schema,
        annoMatrix,
        obsCrossfilter,
      });

      postAsyncSuccessToast("Re-embedding has completed.");
    } catch (error) {
      dispatch({
        type: "reembed: request aborted",
      });
      if (error.name === "AbortError") {
        postAsyncFailureToast("Re-embedding calculation was aborted.");
      } else {
        postNetworkErrorToast(`Re-embedding: ${error.message}`);
      }
      console.log("Reembed exception:", error, error.name, error.message);
    }
  };
}


async function doPreprocessingFetch(dispatch, _getState, reembedParams) {
  const af = abortableFetch(
    `${API.prefix}${API.version}preprocess`,
    {
      method: "PUT",
      headers: new Headers({
        Accept: "application/octet-stream",
        "Content-Type": "application/json",
      }),
      body: JSON.stringify({
        params: reembedParams
      }),
      credentials: "include",
    },
    6000000 // 1 minute timeout
  );
  dispatch({
    type: "preprocess: request start",
    abortableFetch: af,
  });
  const res = await af.ready();

  if (res.ok && res.headers.get("Content-Type").includes("application/json")) {
    // #TOOD: TRIGGER CELL SUBSETTING AND AWAIT RESULTS!
    return res;
  }

  // else an error
  let msg = `Unexpected HTTP response ${res.status}, ${res.statusText}`;
  const body = await res.text();
  if (body && body.length > 0) {
    msg = `${msg} -- ${body}`;
  }
  throw new Error(msg);
}
/*
functions below are dispatch-able
*/
export function requestPreprocessing(reembedParams) {
  return async (dispatch, getState) => {
    try {
      const {
        obsCrossfilter: prevCrossfilter,
        layoutChoice,
      } = getState();
      const res = await doPreprocessingFetch(dispatch, getState, reembedParams);
      const schema = await res.json()

      dispatch({
        type: "preprocess: request completed",
      });

      const baseDataUrl = `${API.prefix}${API.version}`;
      const annoMatrixNew = new AnnoMatrixLoader(baseDataUrl, schema);
      
      const [annoMatrix, obsCrossfilter] = await _switchEmbedding(
        annoMatrixNew,
        prevCrossfilter,
        layoutChoice.current,
        layoutChoice.current
      );
      
      dispatch({
        type: "annoMatrix: init complete",
        annoMatrix,
        obsCrossfilter
      });      
      

      postAsyncSuccessToast("Preprocessing has completed.");
    } catch (error) {
      dispatch({
        type: "preprocess: request aborted",
      });
      if (error.name === "AbortError") {
        postAsyncFailureToast("Preprocessing calculation was aborted.");
      } else {
        postNetworkErrorToast(`Preprocessing: ${error.message}`);
      }
      console.log("Preprocess exception:", error, error.name, error.message);
    }
  };
}
