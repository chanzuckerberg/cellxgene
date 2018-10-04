import _ from "lodash";
import * as kvCache from "../../../src/util/stateManager/keyvalcache";

/*
This is PRIVATE to keyvalcache and must be kept in sync with
any changs ot that module.  Need to Know - to enable error handling test
*/
const cachePrivateKey = "__kvcachekey__";

/*
helper function - promisify setTimeout()
*/
function timeout(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

describe("kvcache API", () => {
  /*
  test the happy path create/set/get API
  */

  test("simple create", () => {
    /* with defaults */
    const kvc = kvCache.create();
    expect(kvc).toBeDefined();
    expect(kvc).toEqual(expect.objectContaining({}));
    expect(kvCache.get(kvc, "test")).toBeUndefined();

    /* with params */
    const kvc1 = kvCache.create(/* lowWatermark */ 99, /* minTTL */ 0);
    expect(kvc1).toBeDefined();
    expect(kvc1).toEqual(expect.objectContaining({}));
  });

  test("set/get", () => {
    /*
    - check basic get/set functionality
    - check set does not mutate source cache
    */
    const keyName = "foo";
    const kvc1 = kvCache.create();
    expect(kvc1).toBeDefined();
    expect(kvCache.get(kvc1, keyName)).toBeUndefined();

    const val2 = [2];
    const kvc2 = kvCache.set(kvc1, keyName, val2);
    expect(kvc2).toBeDefined();
    expect(kvc2).not.toBe(kvc1);
    expect(kvCache.get(kvc1, keyName)).toBeUndefined();
    expect(kvCache.get(kvc2, keyName)).toBe(val2);

    const val3 = [3];
    const kvc3 = kvCache.set(kvc2, keyName, val3);
    expect(kvc3).toBeDefined();
    expect(kvc3).not.toBe(kvc1);
    expect(kvc3).not.toBe(kvc2);
    expect(kvCache.get(kvc1, keyName)).toBeUndefined();
    expect(kvCache.get(kvc2, keyName)).toBe(val2);
    expect(kvCache.get(kvc3, keyName)).toBe(val3);
  });
});

describe("common error handling", () => {
  /*
  Test common error handlers
  */

  test("set() protection from namespace pollution", () => {
    /*
    Test that set() will not allow use of the private cache key
    */
    const kvc = kvCache.create();
    expect(() => {
      kvCache.set(kvc, cachePrivateKey, {});
    }).toThrow();
  });

  test("create() does not accept bogus config", () => {
    expect(() => {
      kvCache.create([], {});
    }).toThrow();
    expect(() => {
      kvCache.create(-99, 0);
    }).toThrow();
    expect(() => {
      kvCache.create(100, -1);
    }).toThrow();
    expect(() => {
      kvCache.create(1000, "foobar");
    }).toThrow();
    expect(() => {
      kvCache.create(null, 8);
    }).toThrow();
  });
});

describe("map", () => {
  /*
  Test kvCache.map() - create new cache that is a transformation of an
  existing cache
  */
  test("map of empty cache", () => {
    const kvc = kvCache.create();
    const callback = jest.fn();
    const kvcMapped = kvCache.map(kvc, callback);
    expect(callback).not.toHaveBeenCalled();
    expect(kvcMapped).toBeDefined();
    expect(kvcMapped).not.toBe(kvc); // immutable operation
    expect(kvcMapped).toEqual(kvc);
  });

  test("map of non-empty cache", () => {
    const key = "aKey";
    const val = [0, 1, 2];
    let kvc = kvCache.create();
    kvc = kvCache.set(kvc, key, val);
    const mockCB = jest.fn().mockImplementation(v => [...v]);
    const kvcMapped = kvCache.map(kvc, mockCB);

    expect(kvcMapped).toBeDefined();
    expect(kvcMapped).not.toBe(kvc); // immutable operation
    expect(_.isEqual(kvc, kvcMapped)).toBe(true);

    expect(mockCB).toHaveBeenCalledTimes(1);
    expect(mockCB).toHaveBeenLastCalledWith(val, key);
  });
});

describe("flush", () => {
  /*
  test various cache flush behavior
  */
  test("flush - lowWatermark, disable minTTL", () => {
    /*
    verify lowWatermark functions correctly
    */

    // set lowWatermark to 2, set three times - only the final two
    // should remain.
    let kvc = kvCache.create(2, 0);
    ["a", "b", "c"].forEach(k => {
      kvc = kvCache.set(kvc, k, []);
    });

    expect(kvc).toEqual(
      expect.objectContaining({
        b: expect.arrayContaining([]),
        c: expect.arrayContaining([])
      })
    );
    expect(kvc).toEqual(
      expect.not.objectContaining({
        a: expect.arrayContaining([])
      })
    );
  });

  test("flush - minTTL, disable lowWatermark", async () => {
    /*
    verify minTTL functions correctly
    */

    // set minTTL to 1 ms
    let kvc = kvCache.create(0, 10);
    kvc = kvCache.set(kvc, "a", []);
    await timeout(20);
    ["b", "c"].forEach(k => {
      kvc = kvCache.set(kvc, k, []);
    });

    expect(kvc).toEqual(
      expect.objectContaining({
        b: expect.arrayContaining([]),
        c: expect.arrayContaining([])
      })
    );
    expect(kvc).toEqual(
      expect.not.objectContaining({
        a: expect.arrayContaining([])
      })
    );
  });

  test("flush - minTTL and lowWatermark", async () => {
    /*
    verify minTTL functions correctly
    */

    // set lowwatermark to 3, minTTL to 1 ms
    let kvc = kvCache.create(3, 10);
    kvc = kvCache.set(kvc, "a", []);
    // delay
    await timeout(20);
    ["b", "c"].forEach(k => {
      kvc = kvCache.set(kvc, k, []);
    });

    expect(kvc).toEqual(
      expect.objectContaining({
        a: expect.arrayContaining([]),
        b: expect.arrayContaining([]),
        c: expect.arrayContaining([])
      })
    );

    kvc = kvCache.set(kvc, "d", []);
    expect(kvc).toEqual(
      expect.objectContaining({
        b: expect.arrayContaining([]),
        c: expect.arrayContaining([]),
        d: expect.arrayContaining([])
      })
    );
    expect(kvc).toEqual(
      expect.not.objectContaining({
        a: expect.arrayContaining([])
      })
    );
  });

  test("manual flush", async () => {
    let kvc = kvCache.create(1, 10);
    ["a", "b", "c", "d"].forEach(k => {
      kvc = kvCache.set(kvc, k, []);
    });

    // Before TTL has expired, should have all values in cache.
    expect(kvc).toEqual(
      expect.objectContaining({
        a: expect.arrayContaining([]),
        b: expect.arrayContaining([]),
        c: expect.arrayContaining([])
      })
    );

    // let TTL expire
    await timeout(10);

    // manually flush
    const postFlushKvc = kvCache.flush(kvc);
    expect(postFlushKvc).toBeDefined();
    expect(postFlushKvc).not.toBe(kvc);
    expect(postFlushKvc).toEqual(
      expect.objectContaining({
        d: expect.arrayContaining([])
      })
    );
  });
});
