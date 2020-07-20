/*

DataframeCache - to improve memoization, we often want to have the same
Dataframe object returned when the contents of the dataframe are identical.
The core dataframe class does not support this, but it does attempt to maintain
strict eq and immutability of columns indices.

This is a function that will maintain a least-recently created cache of Dataframe 
objects.
*/

import { memoize } from "./util";

function hashDataframe(df) {
  if (df.isEmpty()) return "";
  return df.__columnsAccessor.map((c) => c.__id).join(",");
}

function noop(df) {
  return df;
}

const dataframeMemo = (capacity = 100) =>
  memoize(noop, hashDataframe, capacity);

export default dataframeMemo;
