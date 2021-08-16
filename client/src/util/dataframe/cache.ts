/*

DataframeCache - to improve memoization, we often want to have the same
Dataframe object returned when the contents of the dataframe are identical.
The core dataframe class does not support this, but it does attempt to maintain
strict eq and immutability of columns indices.

This is a function that will maintain a least-recently created cache of Dataframe 
objects.
*/

import { memoize } from "./util";
import Dataframe from "./dataframe";

function hashDataframe(df: Dataframe): string {
  if (df.isEmpty()) return "";
  return df.__columnsAccessor.map((c) => c.__id).join(",");
}

function noop(df: Dataframe): Dataframe {
  return df;
}

const dataframeMemo = (capacity = 100): ((df: Dataframe) => Dataframe) =>
  memoize(noop, hashDataframe, capacity);

export default dataframeMemo;
