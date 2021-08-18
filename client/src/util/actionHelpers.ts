import sortBy from "lodash.sortby";
/* XXX: cough, cough, ... */
import { postNetworkErrorToast } from "../components/framework/toasters";
import type { AppDispatch, RootState } from "../reducers";

/*
dispatch an action error to the user.   Currently we use
async toasts.
*/
let networkErrorToastKey: string | null = null;
export const dispatchNetworkErrorMessageToUser = (message: string): void => {
  if (!networkErrorToastKey) {
    networkErrorToastKey = postNetworkErrorToast(message);
  } else {
    postNetworkErrorToast(message, networkErrorToastKey);
  }
};

/*
Catch unexpected errors and make sure we don't lose them!
*/
export function catchErrorsWrap(
  fn: (dispatch: AppDispatch, getState: () => RootState) => Promise<void>,
  dispatchToUser = false
) {
  return (dispatch: AppDispatch, getState: () => RootState): void => {
    fn(dispatch, getState).catch((error: Error) => {
      console.error(error);
      if (dispatchToUser) {
        dispatchNetworkErrorMessageToUser(error.message);
      }
      dispatch({ type: "UNEXPECTED ERROR", error });
    });
  };
}

/**
 * Wrapper to perform async fetch with some modest error handling
 * and decoding.  Arguments are identical to standard fetch.
 */
export const doFetch = async (
  url: string,
  init?: RequestInit
): Promise<Response> => {
  try {
    // add defaults to the fetch init param.
    init = {
      method: "get",
      credentials: "include",
      ...init,
    };
    const acceptType = (init.headers as Headers)?.get("Accept");
    const res = await fetch(url, init);
    if (
      res.ok &&
      (!acceptType || res.headers?.get("Content-Type")?.includes(acceptType))
    ) {
      return res;
    }

    // else an error
    const msg = `Unexpected HTTP response ${res.status}, ${res.statusText}`;
    dispatchNetworkErrorMessageToUser(msg);
    throw new Error(msg);
  } catch (e) {
    // network error
    const msg = "Unexpected HTTP error";
    dispatchNetworkErrorMessageToUser(msg);
    throw e;
  }
};

/*
Wrapper to perform an async fetch and JSON decode response.
*/
export const doJsonRequest = async <T = unknown>(
  url: string,
  init?: RequestInit
): Promise<T> => {
  const res = await doFetch(url, {
    ...init,
    headers: new Headers({ Accept: "application/json" }),
  });
  return res.json();
};

/*
Wrapper to perform an async fetch for binary data.
*/
export const doBinaryRequest = async (
  url: string,
  init?: RequestInit
): Promise<ArrayBuffer> => {
  const res = await doFetch(url, {
    ...init,
    headers: new Headers({ Accept: "application/octet-stream" }),
  });
  return res.arrayBuffer();
};

/*
This function "packs" filter index lists into the more efficient
"range" form specified in the REST 0.2 spec.

Specifically, it turns an array of indices [0, 1, 2, 10, 11, 14, ...]
into a form that encodes runs of consecutive numbers as [min, max].
Array may not be sorted, but will only contain uniq values.

Parameters:
   indices - input array of numbers (index)
   minRangeLength - hint, min range length before it is encoded into range format.
   sorted - boolean hint indicating array is presorted, ascending order

So [1, 2, 3, 4, 10, 11, 14] -> [ [1, 4], [10, 11], 14]
*/
export const rangeEncodeIndices = (
  indices: Array<number>,
  minRangeLength = 3,
  sorted = false
): Array<number | [number, number]> => {
  if (indices.length === 0) {
    return indices;
  }

  if (!sorted) {
    indices = sortBy(indices);
  }

  const result = new Array(indices.length);
  let resultTail = 0;

  let i = 0;
  while (i < indices.length) {
    const begin = indices[i];
    let current;
    do {
      current = indices[i];
      i += 1;
    } while (i < indices.length && indices[i] === current + 1);

    if (current - begin + 1 >= minRangeLength) {
      result[resultTail] = [begin, current];
      resultTail += 1;
    } else {
      for (let j = begin; j <= current; j += 1, resultTail += 1) {
        result[resultTail] = j;
      }
    }
  }

  result.length = resultTail;
  return result;
};
