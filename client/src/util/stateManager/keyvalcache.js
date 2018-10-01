// jshint esversion: 6
import _ from "lodash";

/*
Very simple key/value cache for use by World & Universe.

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

function create(lowWatermark = 32, minTTL = 1000) {
  return {
    [cachePrivateKey]: {
      lowWatermark,
      minTTL
    }
  };
}

function get(kvcache, key) {
  const val = kvcache[key];
  if (val) {
    val[cachePrivateKey] = Date.now();
  }
  return val;
}

function set(kvcache, key, val) {
  const newKvCache = { ...kvcache };
  newKvCache[key] = val;
  val[cachePrivateKey] = Date.now();
  flush(newKvCache, newKvCache[cachePrivateKey].minTTL);
  return newKvCache;
}

/*
Flush elements from cache IF cache size is greater than lowWatermark, and
those elements are older than minAgeMS
*/
function flush(kvcache, minAgeMs = 0) {
  if (minAgeMs < 0) return kvcache;
  const eol = Date.now() - minAgeMs;
  const { lowWatermark } = kvcache[cachePrivateKey];
  const keys = _(kvcache)
    .keys()
    .filter(k => k !== cachePrivateKey)
    .filter(k => kvcache[k][cachePrivateKey] < eol)
    .sortBy([k => kvcache[k][cachePrivateKey]])
    .value();

  if (keys.length > lowWatermark) {
    const numKeysToDelete = keys.length - lowWatermark;
    const keysToDelete = _.slice(keys, 0, numKeysToDelete);
    _.forEach(keysToDelete, k => delete kvcache[k]);
  }

  return kvcache;
}

/*
use to create a cache that is a transformation of another cache.
*/
function map(srcKvCache, cb, createOptions) {
  const keysInSrcKvCache = _(srcKvCache)
    .keys()
    .filter(k => k !== cachePrivateKey)
    .value();
  const newKvCache = create(createOptions.lowWatermark, createOptions.minTTL);
  _.forEach(keysInSrcKvCache, key => {
    const val = cb(get(srcKvCache, key));
    newKvCache[key] = val;
    val[cachePrivateKey] = Date.now();
  });
  return newKvCache;
}

export { create, get, set, flush, map };
