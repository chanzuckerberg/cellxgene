// jshint esversion: 6
import _ from "lodash";

/*
Very simple key/value cache for use by World & Universe.   Cache keys must
be a string, and values are any JS non-primitive value.

  * constructor(lowWatermark, minTTL):
      - lowWatermark defines the number of cache elements below which
        flushing will not occur.
      - minTTL defines minimum time in milliseconds that cache entries will live.
        A value of -1 disables automatic flushing (flush() can still
        be called by external user).
  * set() - add a key/val pair.
  * get() - get a value or undefined if not present.
  * flush(minAgeMs) - flush cache entries in excess of lowWatermark if those
    entries are older than minAgeMs.

*/

const cachePrivateKey = "__kvcachekey__";
const defaultLowWatermark = 32;
const defaultMinTTL = 1000;

function create(lowWatermark = defaultLowWatermark, minTTL = defaultMinTTL) {
  if (typeof minTTL !== "number" || typeof lowWatermark !== "number") {
    throw new TypeError(
      "minTTL and lowWatermark parameters must be a primitive number"
    );
  }
  if (lowWatermark < 0 || minTTL < 0) {
    throw new RangeError(
      "minTTL and lowWatermark parameters must be number greater than zero"
    );
  }

  return {
    [cachePrivateKey]: {
      lowWatermark,
      minTTL
    }
  };
}

function get(kvcache, key) {
  if (key === cachePrivateKey) {
    throw new RangeError(`key parameter may not have value ${cachePrivateKey}`);
  }

  const val = kvcache[key];
  if (val) {
    val[cachePrivateKey] = Date.now();
  }
  return val;
}

function set(kvcache, key, val) {
  if (key === cachePrivateKey) {
    throw new RangeError(`key parameter may not have value ${cachePrivateKey}`);
  }

  const newKvCache = { ...kvcache };
  newKvCache[key] = val;
  val[cachePrivateKey] = Date.now();
  flushInPlace(newKvCache);
  return newKvCache;
}

function flush(kvcache) {
  const newKvCache = { ...kvcache };
  flushInPlace(newKvCache);
  return newKvCache;
}

/*
Flush elements from cache IF cache size is greater than lowWatermark, and
those elements are older than minAgeMS
*/
function flushInPlace(kvCache) {
  const { lowWatermark, minTTL } = kvCache[cachePrivateKey];
  const eol = Date.now() - minTTL;
  const allKeys = _(kvCache)
    .keys()
    .filter(k => k !== cachePrivateKey)
    .sortBy([k => kvCache[k][cachePrivateKey]])
    .value();

  if (allKeys.length > lowWatermark) {
    const keysToDelete = _(allKeys)
      .slice(0, allKeys.length - lowWatermark)
      .filter(k => kvCache[k][cachePrivateKey] <= eol)
      .value();
    _.forEach(keysToDelete, k => delete kvCache[k]);
  }

  return kvCache;
}

/*
use to create a cache that is a transformation of another cache.
*/
function map(srcKvCache, cb, createOptions) {
  const keysInSrcKvCache = _(srcKvCache)
    .keys()
    .filter(k => k !== cachePrivateKey)
    .value();
  const lowWatermark = _.get(
    createOptions,
    "lowWatermark",
    defaultLowWatermark
  );
  const minTTL = _.get(createOptions, "minTTL", defaultMinTTL);
  const newKvCache = create(lowWatermark, minTTL);
  _.forEach(keysInSrcKvCache, key => {
    const val = cb(get(srcKvCache, key), key);
    newKvCache[key] = val;
    val[cachePrivateKey] = Date.now();
  });
  return newKvCache;
}

export { create, get, set, flush, map };
