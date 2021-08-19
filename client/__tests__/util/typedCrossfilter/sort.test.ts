import {
  sortArray,
  sortIndex,
  lowerBound,
} from "../../../src/util/typedCrossfilter/sort";

/*
Sort tests should keep in mind that there are separate code
paths for:
  - small vs. large arrays (insertionsort only)
  - float-only typed arrays vs. other array types (non-finite handling)
  - indexed vs. direct sort
*/

const pInf = Number.POSITIVE_INFINITY;
const nInf = Number.NEGATIVE_INFINITY;

function fillRange(arr: any, start = 0) {
  const larr = arr;
  for (let i = 0, len = larr.length; i < len; i += 1) {
    larr[i] = i + start;
  }
  return larr;
}

function fillRand(arr: any) {
  for (let i = 0, len = arr.length; i < len; i += 1) {
    arr[i] = Math.random();
  }
  return arr;
}

describe("sortArray", () => {
  describe("JS vals", () => {
    [
      [true, false],
      ["a", "b", "0", "1"],
      [0, "a", true, null, undefined, 3.1415],
      fillRand(new Array(1000)),
      ["a", NaN, null, pInf],
    ].map((val, idx) =>
      test(`JS vals ${idx}`, () => {
        expect(sortArray(val)).toMatchObject(val.sort());
      })
    );
  });

  describe("finite numbers", () => {
    [Array, Float32Array, Uint32Array, Int32Array, Float64Array].map((Type) =>
      test(Type.name, () => {
        // @ts-expect-error ts-migrate(2349) FIXME: This expression is not callable.
        expect(sortArray(Type.from([6, 5, 4, 3, 2, 1, 0]))).toMatchObject(
          // @ts-expect-error ts-migrate(2349) FIXME: This expression is not callable.
          Type.from([0, 1, 2, 3, 4, 5, 6])
        );

        // @ts-expect-error ts-migrate(2349) FIXME: This expression is not callable.
        expect(sortArray(Type.from([6, 5, 4, 3, 2, 1]))).toMatchObject(
          // @ts-expect-error ts-migrate(2349) FIXME: This expression is not callable.
          Type.from([1, 2, 3, 4, 5, 6])
        );

        const source = fillRand(new Type(1000));
        // @ts-expect-error ts-migrate(2349) FIXME: This expression is not callable.
        expect(sortArray(Type.from(source))).toMatchObject(
          // @ts-expect-error ts-migrate(2349) FIXME: This expression is not callable.
          Type.from(source).sort()
        );
      })
    );
  });

  describe("non-finite numbers", () => {
    test("infinity", () => {
      expect(sortArray(new Float32Array([pInf, nInf, 0, 1, 2]))).toMatchObject(
        new Float32Array([nInf, 0, 1, 2, pInf])
      );
      expect(
        sortArray(new Float32Array([pInf, nInf, pInf, nInf]))
      ).toMatchObject(new Float32Array([nInf, nInf, pInf, pInf]));
      expect(
        sortArray(new Float32Array([pInf, nInf, pInf, nInf, pInf]))
      ).toMatchObject(new Float32Array([nInf, nInf, pInf, pInf, pInf]));
      expect(
        sortArray(
          new Float32Array(100).fill(Infinity, 0, 50).fill(-Infinity, 50, 100)
        )
      ).toMatchObject(
        new Float32Array(100).fill(-Infinity, 0, 50).fill(Infinity, 50, 100)
      );
    });
    test("NaN", () => {
      expect(sortArray(new Float64Array([NaN, 2, 1, 0]))).toMatchObject(
        new Float64Array([0, 1, 2, NaN])
      );
      expect(sortArray(new Float32Array([NaN, 2, 1, 0]))).toMatchObject(
        new Float32Array([0, 1, 2, NaN])
      );
      expect(sortArray(new Float32Array([NaN, 2, NaN, 1, 0]))).toMatchObject(
        new Float32Array([0, 1, 2, NaN, NaN])
      );
      expect(sortArray(new Float32Array([NaN, 2, 1, NaN, 0]))).toMatchObject(
        new Float32Array([0, 1, 2, NaN, NaN])
      );
      expect(
        sortArray(fillRange(new Float32Array(100)).fill(NaN, 0, 10))
      ).toMatchObject(fillRange(new Float32Array(100), 10).fill(NaN, 90, 100));
    });

    test("mixed numbers", () => {
      expect(
        sortArray(new Float32Array([NaN, pInf, nInf, NaN, NaN]))
      ).toMatchObject(new Float32Array([nInf, pInf, NaN, NaN, NaN]));
      expect(
        sortArray(new Float32Array([NaN, pInf, nInf, NaN, 1, NaN, 2]))
      ).toMatchObject(new Float32Array([nInf, 1, 2, pInf, NaN, NaN, NaN]));
      expect(
        sortArray(new Float32Array([NaN, pInf, nInf, 0, 1, NaN, 2]))
      ).toMatchObject(new Float32Array([nInf, 0, 1, 2, pInf, NaN, NaN]));
      expect(
        sortArray(
          fillRange(new Float32Array(100))
            .fill(NaN, 0, 10)
            .fill(Infinity, 10, 20)
        )
      ).toMatchObject(
        fillRange(new Float32Array(100), 20)
          .fill(Infinity, 80, 90)
          .fill(NaN, 90, 100)
      );
    });
  });
});

describe("sortIndex", () => {
  describe("finite numbers", () => {
    [Array, Float32Array, Uint32Array, Int32Array, Float64Array].map((Type) =>
      test(Type.name, () => {
        // @ts-expect-error ts-migrate(2349) FIXME: This expression is not callable.
        const source1 = Type.from([6, 5, 4, 3, 2, 1, 0]);
        const index1 = fillRange(new Uint32Array(source1.length));
        expect(sortIndex(index1, source1)).toMatchObject(
          index1.sort((a: any, b: any) => source1[a] - source1[b])
        );

        // @ts-expect-error ts-migrate(2349) FIXME: This expression is not callable.
        const source2 = Type.from([6, 5, 4, 3, 2, 1]);
        const index2 = fillRange(new Uint32Array(source2.length));
        expect(sortIndex(index2, source2)).toMatchObject(
          index2.sort((a: any, b: any) => source1[a] - source1[b])
        );

        const source3 = fillRand(new Type(1000));
        const index3 = fillRange(new Uint32Array(source3.length));
        expect(sortIndex(index3, source3)).toMatchObject(
          index3.sort((a: any, b: any) => source1[a] - source1[b])
        );
      })
    );
  });

  describe("non-finite numbers", () => {
    test("mixed numbers", () => {
      const source1 = new Float32Array([NaN, pInf, nInf, NaN, 1, NaN, 2]);
      const index1 = fillRange(new Uint32Array(source1.length));
      expect(sortIndex(index1, source1)).toMatchObject(
        new Uint32Array([2, 4, 6, 1, 0, 3, 5])
      );

      const source2 = new Float32Array([NaN, pInf, nInf, 0, 1, NaN, 2]);
      const index2 = fillRange(new Uint32Array(source2.length));
      expect(sortIndex(index2, source2)).toMatchObject(
        new Uint32Array([2, 3, 4, 6, 1, 0, 5])
      );
    });
  });
});

describe("lowerBound", () => {
  test("non-float path", () => {
    expect(lowerBound([], 0, 0, 0)).toEqual(0);

    expect(lowerBound([0, 1, 2, 3], -1, 0, 4)).toEqual(0);
    expect(lowerBound([0, 1, 2, 3], 0, 0, 4)).toEqual(0);
    expect(lowerBound([0, 1, 2, 3], 1, 0, 4)).toEqual(1);
    expect(lowerBound([0, 1, 2, 3], 3, 0, 4)).toEqual(3);
    expect(lowerBound([0, 1, 2, 3], 4, 0, 4)).toEqual(4);
    expect(lowerBound([0, 1, 2, 3], 4, 0, 3)).toEqual(3);

    expect(lowerBound([0, 1, 2, 3, 4], -1, 0, 5)).toEqual(0);
    expect(lowerBound([0, 1, 2, 3, 4], 0, 0, 5)).toEqual(0);
    expect(lowerBound([0, 1, 2, 3, 4], 2, 0, 5)).toEqual(2);
    expect(lowerBound([0, 1, 2, 3, 4], 4, 0, 5)).toEqual(4);
    expect(lowerBound([0, 1, 2, 3, 4], 5, 0, 5)).toEqual(5);

    expect(lowerBound([0, 2, 4, 6, 8], 5, 0, 5)).toEqual(3);
    expect(lowerBound([0, 2, 2, 2, 8], 5, 0, 5)).toEqual(4);

    expect(lowerBound([0, 1, 2, 3, 4, 5, 6, 7, 8], 3, 2, 4)).toEqual(3);
    expect(lowerBound([0, 1, 2, 3, 4, 5, 6, 7, 8], 99, 2, 4)).toEqual(4);
  });

  test("float path, finites", () => {
    expect(lowerBound(new Float32Array([0, 1, 2, 3]), 1, 0, 4)).toEqual(1);

    expect(lowerBound(new Float32Array([]), 0, 0, 0)).toEqual(0);

    expect(lowerBound(new Float32Array([0, 1, 2, 3]), -1, 0, 4)).toEqual(0);
    expect(lowerBound(new Float32Array([0, 1, 2, 3]), 0, 0, 4)).toEqual(0);
    expect(lowerBound(new Float32Array([0, 1, 2, 3]), 1, 0, 4)).toEqual(1);
    expect(lowerBound(new Float32Array([0, 1, 2, 3]), 3, 0, 4)).toEqual(3);
    expect(lowerBound(new Float32Array([0, 1, 2, 3]), 4, 0, 4)).toEqual(4);
    expect(lowerBound(new Float32Array([0, 1, 2, 3]), 4, 0, 3)).toEqual(3);

    expect(lowerBound(new Float32Array([0, 1, 2, 3, 4]), -1, 0, 5)).toEqual(0);
    expect(lowerBound(new Float32Array([0, 1, 2, 3, 4]), 0, 0, 5)).toEqual(0);
    expect(lowerBound(new Float32Array([0, 1, 2, 3, 4]), 2, 0, 5)).toEqual(2);
    expect(lowerBound(new Float32Array([0, 1, 2, 3, 4]), 4, 0, 5)).toEqual(4);
    expect(lowerBound(new Float32Array([0, 1, 2, 3, 4]), 5, 0, 5)).toEqual(5);

    expect(lowerBound(new Float32Array([0, 2, 4, 6, 8]), 5, 0, 5)).toEqual(3);
    expect(lowerBound(new Float32Array([0, 2, 2, 2, 8]), 5, 0, 5)).toEqual(4);

    expect(
      lowerBound(new Float32Array([0, 1, 2, 3, 4, 5, 6, 7, 8]), 3, 2, 4)
    ).toEqual(3);
    expect(
      lowerBound(new Float32Array([0, 1, 2, 3, 4, 5, 6, 7, 8]), 99, 2, 4)
    ).toEqual(4);
  });

  test("float path, non-finite", () => {
    expect(
      lowerBound(
        new Float32Array([-Infinity, 0, 1, Infinity, NaN]),
        -Infinity,
        0,
        5
      )
    ).toEqual(0);
    expect(
      lowerBound(new Float32Array([-Infinity, 0, 1, Infinity, NaN]), 0, 0, 5)
    ).toEqual(1);
    expect(
      lowerBound(new Float32Array([-Infinity, 0, 1, Infinity, NaN]), 1, 0, 5)
    ).toEqual(2);
    expect(
      lowerBound(new Float32Array([-Infinity, 0, 1, Infinity, NaN]), 2, 0, 5)
    ).toEqual(3);
    expect(
      lowerBound(
        new Float32Array([-Infinity, 0, 1, Infinity, NaN]),
        Infinity,
        0,
        5
      )
    ).toEqual(3);
    expect(
      lowerBound(new Float32Array([-Infinity, 0, 1, Infinity, NaN]), NaN, 0, 5)
    ).toEqual(4);
  });
});
