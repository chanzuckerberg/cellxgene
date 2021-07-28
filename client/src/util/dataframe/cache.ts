/*

DataframeCache - to improve memoization, we often want to have the same
Dataframe object returned when the contents of the dataframe are identical.
The core dataframe class does not support this, but it does attempt to maintain
strict eq and immutability of columns indices.

This is a function that will maintain a least-recently created cache of Dataframe 
objects.
*/

import { memoize } from "./util";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function hashDataframe(df: any) {
  if (df.isEmpty()) return "";
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  return df.__columnsAccessor.map((c: any) => c.__id).join(",");
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function noop(df: any) {
  return df;
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
const dataframeMemo = (capacity = 100) =>
  memoize(noop, hashDataframe, capacity);

export default dataframeMemo;
