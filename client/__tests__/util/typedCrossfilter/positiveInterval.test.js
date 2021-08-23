// const PositiveIntervals = require("../../src/util/typedCrossfilter/positiveIntervals");
import PositiveIntervals from "../../../src/util/typedCrossfilter/positiveIntervals";

describe("canonicalize", () => {
  test("empty", () => {
    expect(PositiveIntervals.canonicalize([])).toEqual([]);
  });

  test("simple, already correct", () => {
    expect(PositiveIntervals.canonicalize([[0, 1]])).toEqual([[0, 1]]);
    expect(
      PositiveIntervals.canonicalize([
        [0, 1],
        [2, 3],
      ])
    ).toEqual([
      [0, 1],
      [2, 3],
    ]);
  });

  test("non-canonical, need to be canonicalized", () => {
    expect(
      PositiveIntervals.canonicalize([
        [0, 1],
        [1, 2],
      ])
    ).toEqual([[0, 2]]);
    expect(
      PositiveIntervals.canonicalize([
        [1, 2],
        [2, 3],
      ])
    ).toEqual([[1, 3]]);
  });
});

describe("union", () => {
  test("empty range", () => {
    expect(PositiveIntervals.union([], [])).toEqual([]);
    expect(PositiveIntervals.union([], [[1, 2]])).toEqual([[1, 2]]);
    expect(
      PositiveIntervals.union(
        [],
        [
          [1, 2],
          [3, 4],
        ]
      )
    ).toEqual([
      [1, 2],
      [3, 4],
    ]);
    expect(PositiveIntervals.union([[3, 4]], [])).toEqual([[3, 4]]);
    expect(
      PositiveIntervals.union(
        [
          [1, 2],
          [3, 4],
        ],
        []
      )
    ).toEqual([
      [1, 2],
      [3, 4],
    ]);
    expect(PositiveIntervals.union([[3, 3]], [])).toEqual([[3, 3]]);
    expect(PositiveIntervals.union([], [[3, 3]])).toEqual([[3, 3]]);
  });

  test("simple ranges", () => {
    expect(PositiveIntervals.union([[1, 2]], [[2, 3]])).toEqual([[1, 3]]);
    expect(PositiveIntervals.union([[2, 3]], [[1, 2]])).toEqual([[1, 3]]);
    expect(PositiveIntervals.union([[1, 2]], [[3, 4]])).toEqual([
      [1, 2],
      [3, 4],
    ]);
    expect(
      PositiveIntervals.union(
        [
          [1, 2],
          [3, 4],
        ],
        [
          [6, 7],
          [19, 40],
        ]
      )
    ).toEqual([
      [1, 2],
      [3, 4],
      [6, 7],
      [19, 40],
    ]);
    expect(
      PositiveIntervals.union(
        [[1, 4]],
        [
          [1, 1],
          [3, 4],
        ]
      )
    ).toEqual([[1, 4]]);
    expect(PositiveIntervals.union([[3, 3]], [[4, 4]])).toEqual([
      [3, 3],
      [4, 4],
    ]);
  });
});

describe("intersection", () => {
  test("empty range", () => {
    expect(PositiveIntervals.intersection([], [])).toEqual([]);
    expect(PositiveIntervals.intersection([], [[1, 2]])).toEqual([]);
    expect(PositiveIntervals.intersection([[1, 2]], [])).toEqual([]);
  });

  test("simple", () => {
    expect(PositiveIntervals.intersection([[1, 2]], [[2, 3]])).toEqual([]);
    expect(PositiveIntervals.intersection([[2, 3]], [[1, 2]])).toEqual([]);
    expect(PositiveIntervals.intersection([[1, 10]], [[1, 10]])).toEqual([
      [1, 10],
    ]);
    expect(PositiveIntervals.intersection([[1, 10]], [[2, 8]])).toEqual([
      [2, 8],
    ]);
    expect(PositiveIntervals.intersection([[2, 8]], [[1, 10]])).toEqual([
      [2, 8],
    ]);
    expect(PositiveIntervals.intersection([[1, 10]], [[2, 12]])).toEqual([
      [2, 10],
    ]);
    expect(PositiveIntervals.intersection([[2, 12]], [[1, 10]])).toEqual([
      [2, 10],
    ]);
    expect(PositiveIntervals.intersection([[1, 10]], [[1, 8]])).toEqual([
      [1, 8],
    ]);
    expect(PositiveIntervals.intersection([[1, 8]], [[1, 10]])).toEqual([
      [1, 8],
    ]);
    expect(
      PositiveIntervals.intersection(
        [[1, 10]],
        [
          [1, 2],
          [6, 9],
        ]
      )
    ).toEqual([
      [1, 2],
      [6, 9],
    ]);
    expect(PositiveIntervals.intersection([[0, 2638]], [[1363, 2638]])).toEqual(
      [[1363, 2638]]
    );
    expect(PositiveIntervals.intersection([[1, 2]], [[1, 2]])).toEqual([
      [1, 2],
    ]);
  });
});

describe("difference", () => {
  test("empty", () => {
    expect(PositiveIntervals.difference([], [])).toEqual([]);
    expect(PositiveIntervals.difference([], [[1, 10]])).toEqual([]);
    expect(PositiveIntervals.difference([[1, 10]], [])).toEqual([[1, 10]]);
  });

  test("simple", () => {
    expect(
      PositiveIntervals.difference(
        [
          [1, 2],
          [3, 4],
        ],
        []
      )
    ).toEqual([
      [1, 2],
      [3, 4],
    ]);
    expect(
      PositiveIntervals.difference(
        [
          [1, 2],
          [3, 10],
        ],
        [[5, 10]]
      )
    ).toEqual([
      [1, 2],
      [3, 5],
    ]);
    expect(
      PositiveIntervals.difference(
        [
          [1, 2],
          [3, 10],
        ],
        [[0, 5]]
      )
    ).toEqual([[5, 10]]);
    expect(
      PositiveIntervals.difference(
        [[0, 2638]],
        [
          [0, 1363],
          [2055, 2638],
        ]
      )
    ).toEqual([[1363, 2055]]);
    expect(
      PositiveIntervals.difference(
        [
          [0, 1363],
          [2055, 2638],
        ],
        [[0, 2638]]
      )
    ).toEqual([]);
    expect(PositiveIntervals.difference([[0, 10]], [[0, 1]])).toEqual([
      [1, 10],
    ]);
    expect(PositiveIntervals.difference([[0, 10]], [[1, 2]])).toEqual([
      [0, 1],
      [2, 10],
    ]);
    expect(PositiveIntervals.difference([[0, 10]], [[9, 10]])).toEqual([
      [0, 9],
    ]);
  });
});
