import * as kvCache from "../../../src/util/stateManager/keyvalcache";

describe("kvcache API", () => {
  test("simple create", () => {
    const kvc = kvCache.create();
    expect(kvc).toBeDefined();
    expect(kvc).toEqual(expect.objectContaining({}));
    expect(kvCache.get(kvc, "test")).toBeUndefined();
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

/*
tests to write:

- map
- flush
- flush on set if limits are hit

*/
