import _ from "lodash";
/* XXX: cough, cough, ... */
import { postNetworkErrorToast } from "../components/framework/toasters";

/*
dispatch an action error to the user.   Currently we use
async toasts.
*/
export const dispatchNetworkErrorMessageToUser = message =>
  postNetworkErrorToast(message);

/*
Catch unexpected errors and make sure we don't lose them!
*/
export function catchErrorsWrap(fn, dispatchToUser = false) {
  return (dispatch, getState) => {
    fn(dispatch, getState).catch(error => {
      console.error(error);
      if (dispatchToUser) {
        dispatchNetworkErrorMessageToUser(error.message);
      }
      dispatch({ type: "UNEXPECTED ERROR", error });
    });
  };
}

/*
Wrapper to perform async fetch with some modest error handling
and decoding.
*/
const doFetch = async (url, acceptType) => {
  const res = await fetch(url, {
    method: "get",
    headers: new Headers({
      Accept: acceptType
    })
  });
  if (res.ok && res.headers.get("Content-Type").includes(acceptType)) {
    return res;
  }
  // else an error
  let msg = `Unexpected HTTP response ${res.status}, ${res.statusText}`;
  const body = await res.text();
  if (body && body.length > 0) {
    msg = `${msg} -- ${body}`;
  }
  dispatchNetworkErrorMessageToUser(msg);
  throw new Error(msg);
};

/*
Wrapper to perform an async fetch and JSON decode response.
*/
export const doJsonRequest = async url => {
  const res = await doFetch(url, "application/json");
  return res.json();
};

/*
Wrapper to perform an async fetch for binary data.
*/
export const doBinaryRequest = async url => {
  const res = await doFetch(url, "application/octet-stream");
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
  indices,
  minRangeLength = 3,
  sorted = false
) => {
  if (indices.length === 0) {
    return indices;
  }

  if (!sorted) {
    indices = _.sortBy(indices);
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
