import { API } from "../globals";
import {
  postNetworkErrorToast,
  postAsyncSuccessToast,
  postAsyncFailureToast,
} from "../components/framework/toasters";


async function doSankeyFetch(dispatch, getState, sankeyCategories) {
    const state = getState();
  
    const af = abortableFetch(
      `${API.prefix}${API.version}layout/obs`,
      {
        method: "PUT",
        headers: new Headers({
          Accept: "application/octet-stream",
          "Content-Type": "application/json",
        }),
        body: JSON.stringify({
          method: "sankey",
          sankeyCategories: sankeyCategories,
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

export function requestSankey(sankeyCategories) {
    return async (dispatch, getState) => {
      try {
        const res = await doSankeyFetch(dispatch, getState, sankeyCategories);
        const sankey = await res.json();
        dispatch({
          type: "sankey: request completed",
        });

        postAsyncSuccessToast("Sankey has completed.");

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
  