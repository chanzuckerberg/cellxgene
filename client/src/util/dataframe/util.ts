/*
Private utility code for dataframe
*/
import { AnyFunction, HashArgsFunction } from "./types";

export function callOnceLazy<T>(f: AnyFunction<T>): AnyFunction<T> {
  /*
  call function once, and save the result, regardless of arguments (this is not
  the same as typical memoization).
  */
  let value: T;
  let calledOnce = false;
  const result = function result(...args: unknown[]): T {
    if (!calledOnce) {
      value = f(...args);
      calledOnce = true;
    }
    return value;
  };
  return result;
}

export function memoize<T>(
  fn: AnyFunction<T>,
  hashFn: HashArgsFunction,
  maxResultsCached = -1
): AnyFunction<T> {
  /* 
  function memoization, with user-provided hash.  hashFn must return a
  key which will be unique as a Map key (ie, obeys "sameValueZero" algorithm
  as defined in the JS spec).  For more info on hash key, see:
  https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Map#Key_equality
  */
  const cache = new Map();
  const wrap = function wrap(...args: unknown[]) {
    const key = hashFn(...args);
    if (cache.has(key)) {
      return cache.get(key);
    }
    const result = fn(...args);
    cache.set(key, result);

    if (maxResultsCached > -1 && cache.size > maxResultsCached) {
      /* Least recent insertion deletion */
      cache.delete(cache.keys().next().value);
    }

    return result;
  };

  wrap.clear = function clear() {
    /* clear memoization cache */
    cache.clear();
  };

  return wrap;
}

/**
 *memoization helpers - just a global counter.
 */
let __DataframeMemoId__ = 0;
export function __getMemoId(): string {
  const id = __DataframeMemoId__;
  __DataframeMemoId__ += 1;
  return id.toString();
}
