import quantile from "../../src/util/quantile";

describe("quantile", () => {
  test("single q", () => {
    const arr = new Float32Array([9, 3, 5, 6, 0]);
    expect(quantile([1.0], arr)).toMatchObject([9]);
    expect(quantile([0.9], arr)).toMatchObject([9]);
    expect(quantile([0.8], arr)).toMatchObject([9]);
    expect(quantile([0.7], arr)).toMatchObject([6]);
    expect(quantile([0.6], arr)).toMatchObject([6]);
    expect(quantile([0.5], arr)).toMatchObject([5]);
    expect(quantile([0.4], arr)).toMatchObject([5]);
    expect(quantile([0.3], arr)).toMatchObject([3]);
    expect(quantile([0.2], arr)).toMatchObject([3]);
    expect(quantile([0.1], arr)).toMatchObject([0]);
    expect(quantile([0], arr)).toMatchObject([0]);
  });

  test("multi q", () => {
    const arr = new Float32Array([9, 3, 5, 6, 0]);
    expect(quantile([0, 0.25, 0.5, 0.75, 1.0], arr)).toMatchObject([
      0, 3, 5, 6, 9,
    ]);
  });
});
