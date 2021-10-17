import { API } from "../globals";
import {
  postNetworkErrorToast,
  postAsyncFailureToast,
} from "../components/framework/toasters";

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

async function doLeidenFetch(dispatch, getState) {
    const state = getState();
    // get current embedding
    const { layoutChoice, Leiden } = state;
    const { res } = Leiden;
    const name = `leiden_${layoutChoice.current}_res${Math.round((res+Number.EPSILON)*1000)/1000.0}`
    const af = abortableFetch(
      `${API.prefix}${API.version}leiden`,
      {
        method: "PUT",
        headers: new Headers({
          Accept: "application/octet-stream",
          "Content-Type": "application/json",
        }),
        body: JSON.stringify({
          name: layoutChoice.current,
          cName: name,
          resolution: res,
        }),
        credentials: "include",
      },
      600000 // 1 minute timeout
    );
    dispatch({
      type: "leiden: request start",
      abortableFetch: af,
    });
    const result = await af.ready();
  
    if (result.ok && result.headers.get("Content-Type").includes("application/json")) {
      return result;
    }
  
    // else an error
    let msg = `Unexpected HTTP response ${result.status}, ${result.statusText}`;
    const body = await result.text();
    if (body && body.length > 0) {
      msg = `${msg} -- ${body}`;
    }
    throw new Error(msg);
}

export function requestLeiden() {
    return async (dispatch, getState) => {
      try {
        const res = await doLeidenFetch(dispatch, getState);
        const leiden = await res.json();
        dispatch({
          type: "leiden: request completed",
        });

        return leiden      
      } catch (error) {
        dispatch({
          type: "leiden: request aborted",
        });
        if (error.name === "AbortError") {
          postAsyncFailureToast("Leiden clustering was aborted.");
        } else {
          postNetworkErrorToast(`Leiden: ${error.message}`);
        }
        console.log("Leiden exception:", error, error.name, error.message);
      }
    };
}
  