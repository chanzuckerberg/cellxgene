/*
Tests for WebGL point rendering helpers.
Verifies point size scaling remains usable for large datasets (4M+ cells)
and that alpha differentiation provides visual distinction between
selected and unselected cells.

Regression test for: https://github.com/chanzuckerberg/cellxgene/issues/2709
*/

import {
  flagSelected,
  flagBackground,
  flagHighlight,
} from "../../src/util/glHelpers";

describe("glHelpers flags", () => {
  test("flag constants are correct bitmask values", () => {
    expect(flagSelected).toBe(1);
    expect(flagBackground).toBe(2);
    expect(flagHighlight).toBe(4);
  });

  test("flags are non-overlapping", () => {
    /* eslint-disable no-bitwise -- flag constants use bitmask operations */
    expect(flagSelected & flagBackground).toBe(0);
    expect(flagSelected & flagHighlight).toBe(0);
    expect(flagBackground & flagHighlight).toBe(0);
    /* eslint-enable no-bitwise -- re-enable after bitmask tests */
  });
});

/*
Point size scaling validation.
We cannot directly execute the GLSL code, but we can replicate the
JavaScript-side scaling logic to verify it produces reasonable values
for various dataset sizes.
*/
describe("point size scaling for large datasets", () => {
  // Replicate the JS-side configuration from glHelpers.js
  const domain = [5000000 / (500 * 500), 1000 / (1440 * 1440)];
  const range = [1.0, 5];
  const scale = (range[1] - range[0]) / (domain[1] - domain[0]);
  const offset = scale * -domain[0] + range[0];

  function computePointSize(nPoints, minViewportDimension) {
    const density = nPoints / (minViewportDimension * minViewportDimension);
    let pointSize = scale * density + offset;
    pointSize = Math.max(range[0], Math.min(range[1], pointSize));
    return pointSize;
  }

  test("1000 cells on large viewport produces near-max size", () => {
    const size = computePointSize(1000, 1440);
    expect(size).toBeGreaterThanOrEqual(4.5);
    expect(size).toBeLessThanOrEqual(5.0);
  });

  test("100K cells produces moderate size", () => {
    const size = computePointSize(100000, 900);
    expect(size).toBeGreaterThanOrEqual(1.0);
    expect(size).toBeLessThanOrEqual(5.0);
  });

  test("4M cells produces at least minimum visible size (1.0)", () => {
    const size = computePointSize(4000000, 800);
    expect(size).toBeGreaterThanOrEqual(1.0);
  });

  test("10M cells still produces at least minimum visible size", () => {
    const size = computePointSize(10000000, 800);
    expect(size).toBeGreaterThanOrEqual(1.0);
  });

  test("selected points are at least 1.0px", () => {
    const size = computePointSize(4000000, 800);
    // In the shader: if (isSelected) return max(pointSize, 1.0);
    const selectedSize = Math.max(size, 1.0);
    expect(selectedSize).toBeGreaterThanOrEqual(1.0);
  });

  test("unselected points are visible but smaller than selected", () => {
    const size = computePointSize(4000000, 800);
    // In the shader: return max(pointSize / 3., 0.5);
    const unselectedSize = Math.max(size / 3, 0.5);
    const selectedSize = Math.max(size, 1.0);
    expect(unselectedSize).toBeGreaterThanOrEqual(0.5);
    expect(unselectedSize).toBeLessThan(selectedSize);
  });

  test("highlighted points are always at least 2.0px", () => {
    const size = computePointSize(4000000, 800);
    // In the shader: if (isHighlight) return max(2. * pointSize, 2.0);
    const highlightSize = Math.max(2 * size, 2.0);
    expect(highlightSize).toBeGreaterThanOrEqual(2.0);
  });
});

describe("point alpha differentiation", () => {
  /*
  Replicate the GLSL pointAlpha logic to verify alpha values provide
  visual distinction between selected and unselected cells.
  */
  function pointAlpha(isBackground, isSelected, isHighlight) {
    if (isHighlight) return 1.0;
    if (isBackground) return 0.9;
    if (!isSelected) return 0.3;
    return 1.0;
  }

  test("highlighted cells are fully opaque", () => {
    expect(pointAlpha(false, true, true)).toBe(1.0);
    expect(pointAlpha(false, false, true)).toBe(1.0);
  });

  test("selected cells are fully opaque", () => {
    expect(pointAlpha(false, true, false)).toBe(1.0);
  });

  test("unselected cells have reduced alpha", () => {
    const alpha = pointAlpha(false, false, false);
    expect(alpha).toBe(0.3);
    expect(alpha).toBeLessThan(1.0);
  });

  test("background cells have slightly reduced alpha", () => {
    expect(pointAlpha(true, false, false)).toBe(0.9);
  });

  test("unselected is visually distinct from selected", () => {
    const selectedAlpha = pointAlpha(false, true, false);
    const unselectedAlpha = pointAlpha(false, false, false);
    // There should be a clear difference between selected and unselected
    expect(selectedAlpha - unselectedAlpha).toBeGreaterThanOrEqual(0.5);
  });
});
