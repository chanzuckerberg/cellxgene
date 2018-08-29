// jshint esversion: 6
import _ from "lodash";

/*
Very simple key/value cache for use by World & Universe.

  * constructor(lowWatermark, cachekey):
      - lowWatermark defines the number of cache elements below which
        flushing will not occur.
      - minTTL defines minimum time in MS that cache entries will live.
        A value of -1 disables automatic flushing (flush() can still
        be called by external user).
      - cachekey is a key that will be assigned to any value to track age
  * set() - add a key/val pair.
  * get() - get a value or undefined if not present.
  * flush(minAgeMs) - flush cache entries in excess of lowWatermark if those
    entries are older than minAgeMs.

*/

class KeyValCache {
  constructor(lowWatermark = 32, minTTL = -1, cachekey = "__cachekey__") {
    this.cachekey = cachekey;
    this.cache = {};

    /* don't flush if cache size below this number */
    this.lowWatermark = lowWatermark;
    this.minTTL = minTTL;
  }

  set(key, val) {
    this.cache[key] = val;
    val[this.cachekey] = Date.now();
    this.flush(this.minTTL);
    return val;
  }

  get(key) {
    const val = this.cache[key];
    if (val) {
      val[this.cachekey] = Date.now();
    }
    return val;
  }

  /*
  Flush elements from cache IF cache size is greater than lowWatermark, and
  those elements are older than minAgeMS
  */
  flush(minAgeMs = 0) {
    if (minAgeMs < 0) return this;
    const eol = Date.now() - minAgeMs;
    const { cache, cachekey } = this;
    const keys = _(this.cache)
      .keys()
      .filter(k => cache[k][cachekey] < eol)
      .sortBy([k => cache[k][cachekey]])
      .value();

    if (keys.length > this.lowWatermark) {
      const numKeysToDelete = keys.length - this.lowWatermark;
      const keysToDelete = _.slice(keys, 0, numKeysToDelete);
      _.forEach(keysToDelete, k => delete this.cache[k]);
    }

    return this;
  }
}

export default KeyValCache;
