import { API } from "../globals";
import {
  postNetworkErrorToast,
  postAsyncSuccessToast,
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

async function doSankeyFetch(dispatch, getState) {
    const state = getState();
    // get current embedding
    const { layoutChoice, sankeySelection, annoMatrix } = state;
    const { categories } = sankeySelection;
    const labels = []
    for (const [key, value] of Object.entries(categories)) {
      if(value){
        let t = await annoMatrix.fetch("obs",key)
        labels.push(t)
      }
    }
    if (labels.length === 1){
      labels.push(labels[0])
    }

    const af = abortableFetch(
      `${API.prefix}${API.version}sankey`,
      {
        method: "PUT",
        headers: new Headers({
          Accept: "application/octet-stream",
          "Content-Type": "application/json",
        }),
        body: JSON.stringify({
          name: layoutChoice.current,
          labels: labels,
        }),
        credentials: "include",
      },
      60000 // 1 minute timeout
    );
    dispatch({
      type: "sankey: request start",
      abortableFetch: af,
    });
    const res = await af.ready();
  
    if (res.ok && res.headers.get("Content-Type").includes("application/json")) {
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

export function requestSankey() {
    return async (dispatch, getState) => {
      try {

        const res = await doSankeyFetch(dispatch, getState);
        const sankey = await res.json();
        dispatch({
          type: "sankey: request completed",
        });

        return sankey      
      } catch (error) {
        dispatch({
          type: "sankey: request aborted",
        });
        if (error.name === "AbortError") {
          postAsyncFailureToast("Sankey calculation was aborted.");
        } else {
          postNetworkErrorToast(`Sankey: ${error.message}`);
        }
        console.log("Sankey exception:", error, error.name, error.message);
      }
    };
}
  