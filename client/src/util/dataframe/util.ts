/*
Private utility code for dataframe
*/

export { isTypedArray, isArrayOrTypedArray } from "../typeHelpers";

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function callOnceLazy(f: any) {
  /*
  call function once, and save the result, regardless of arguments (this is not
  the same as typical memoization).
  */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  let value: any;
  let calledOnce = false;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const result = function result(...args: any[]) {
    if (!calledOnce) {
      value = f(...args);
      calledOnce = true;
    }
    return value;
  };
  return result;
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function memoize(fn: any, hashFn: any, maxResultsCached = -1) {
  /* 
  function memoization, with user-provided hash.  hashFn must return a
  key which will be unique as a Map key (ie, obeys "sameValueZero" algorithm
  as defined in the JS spec).  For more info on hash key, see:
  https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Map#Key_equality
  */
  const cache = new Map();
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const wrap = function wrap(...args: any[]) {
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
memoization helpers - just a global counter.
**/
let __DataframeMemoId__ = 0;
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function __getMemoId() {
  const id = __DataframeMemoId__;
  __DataframeMemoId__ += 1;
  return id;
}
