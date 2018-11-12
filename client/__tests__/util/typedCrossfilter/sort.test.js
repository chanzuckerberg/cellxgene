import { sort, sortIndex } from "../../../src/util/typedCrossfilter/sort";

function fillRange(arr, start = 0) {
  const larr = arr;
  for (let i = 0, len = larr.length; i < len; i += 1) {
    larr[i] = i + start;
  }
  return larr;
}

function fillRand(arr) {
  for (let i = 0, len = arr.length; i < len; i += 1) {
    arr[i] = Math.random();
  }
  return arr;
}

describe("sort", () => {
  [Array, Float32Array, Uint32Array, Int32Array, Float64Array].map(Type =>
    test(Type.name, () => {
      expect(sort(Type.from([6, 5, 4, 3, 2, 1, 0]))).toMatchObject(
        Type.from([0, 1, 2, 3, 4, 5, 6])
      );
      expect(sort(Type.from([6, 5, 4, 3, 2, 1]))).toMatchObject(
        Type.from([1, 2, 3, 4, 5, 6])
      );

      const source = fillRand(new Type(1000));
      expect(sort(Type.from(source))).toMatchObject(Type.from(source).sort());
    })
  );
});

describe("sortIndex", () => {
  [Array, Float32Array, Uint32Array, Int32Array, Float64Array].map(Type =>
    test(Type.name, () => {
      const source1 = Type.from([6, 5, 4, 3, 2, 1, 0]);
      const index1 = fillRange(new Uint32Array(source1.length));
      expect(sortIndex(index1, source1)).toMatchObject(
        index1.sort((a, b) => source1[a] - source1[b])
      );

      const source2 = Type.from([6, 5, 4, 3, 2, 1]);
      const index2 = fillRange(new Uint32Array(source2.length));
      expect(sortIndex(index2, source2)).toMatchObject(
        index2.sort((a, b) => source1[a] - source1[b])
      );

      const source3 = fillRand(new Type(1000));
      const index3 = fillRange(new Uint32Array(source3.length));
      expect(sortIndex(index3, source3)).toMatchObject(
        index3.sort((a, b) => source1[a] - source1[b])
      );
    })
  );
});
