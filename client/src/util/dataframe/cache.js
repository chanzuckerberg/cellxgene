/*

DataframeCache - to improve memoization, we often want to have the same
Dataframe object returned when the contents of the dataframe are identical.
The core dataframe class does not support this, but it does attempt to maintain
strict eq and immutability of columns indices.

This is a class that will maintain an LRU cache of Dataframe objects.

The only public method, get(), will return a dataframe with strict eq on 
colummns, indices and dimensions.
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

