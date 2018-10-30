import _ from "lodash";
/* XXX: cough, cough, ... */
import { ErrorToastTopCenter } from "../components/framework/toasters";

/*
dispatch an action error to the user.   Currently we use
async toasts.
*/
export const dispatchErrorMessageToUser = message =>
  ErrorToastTopCenter.show({ message });

/*
Catch unexpected errors and make sure we don't lose them!
*/
export function catchErrorsWrap(fn, dispatchToUser = false) {
  return (dispatch, getState) => {
    fn(dispatch, getState).catch(error => {
      console.error(error);
      if (dispatchToUser) {
        dispatchErrorMessageToUser(error.message);
      }
      dispatch({ type: "UNEXPECTED ERROR", error });
    });
  };
}

/*
Wrapper to perform an async fetch and JSON decode response.
*/
export const doJsonRequest = async url => {
  const res = await fetch(url, {
    method: "get",
    headers: new Headers({
      "Content-Type": "application/json",
      "Accept-Encoding": "gzip, deflate, br"
    })
  });
  return res.json();
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
