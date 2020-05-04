import { rangeEncodeIndices } from "../../src/util/actionHelpers";

describe("rangeEncodeIndices", () => {
  test("small array edge cases", () => {
    expect(rangeEncodeIndices([])).toMatchObject([]);
    expect(rangeEncodeIndices([1])).toMatchObject([1]);
    expect(rangeEncodeIndices([1, 99])).toMatchObject([1, 99]);
    expect(rangeEncodeIndices([99, 1])).toMatchObject([1, 99]);
    expect(rangeEncodeIndices([1, 100, 4])).toMatchObject([1, 4, 100]);
  });

  test("sorted flag", () => {
    expect(rangeEncodeIndices([1, 9, 432], 10, true)).toMatchObject([
      1,
      9,
      432,
    ]);
    expect(rangeEncodeIndices([1, 9, 432], 10, false)).toMatchObject([
      1,
      9,
      432,
    ]);
    expect(
      rangeEncodeIndices([0, 1, 2, 3, 9, 10, 432], 2, true)
    ).toMatchObject([[0, 3], [9, 10], 432]);
    expect(
      rangeEncodeIndices([0, 1, 2, 3, 9, 10, 432], 2, false)
    ).toMatchObject([[0, 3], [9, 10], 432]);
  });

  test("begin or end edge cases", () => {
    expect(
      rangeEncodeIndices([3, 4, 5, 6, 7, 10, 11, 12, 99], 3, false)
    ).toMatchObject([[3, 7], [10, 12], 99]);
    expect(
      rangeEncodeIndices([3, 4, 5, 6, 7, 10, 11, 12], 3, false)
    ).toMatchObject([
      [3, 7],
      [10, 12],
    ]);
    expect(
      rangeEncodeIndices([0, 3, 4, 5, 6, 7, 10, 11, 12], 3, false)
    ).toMatchObject([0, [3, 7], [10, 12]]);
    expect(
      rangeEncodeIndices([0, 3, 4, 5, 6, 7, 10, 11, 12, 99], 3, false)
    ).toMatchObject([0, [3, 7], [10, 12], 99]);
  });

  test("minRangeLength", () => {
    expect(
      rangeEncodeIndices([3, 4, 5, 6, 7, 10, 11, 12, 99], 4, false)
    ).toMatchObject([[3, 7], 10, 11, 12, 99]);
    expect(
      rangeEncodeIndices([3, 4, 5, 6, 7, 10, 11, 12], 4, false)
    ).toMatchObject([[3, 7], 10, 11, 12]);
    expect(
      rangeEncodeIndices([0, 3, 4, 5, 6, 7, 10, 11, 12], 4, false)
    ).toMatchObject([0, [3, 7], 10, 11, 12]);
    expect(
      rangeEncodeIndices([0, 3, 4, 5, 6, 7, 10, 11, 12, 99], 4, false)
    ).toMatchObject([0, [3, 7], 10, 11, 12, 99]);
  });
});
