/*
Utility code for WebGL shaders
*/

/*
PointFlags:
  Point flags are used in graph & scatter plots.

  We want a bitmask-like flag structure, but due to webgl limitations
  must emulate it with floats.

  Supported flags are:

    selected: the point is currently selected
    highlight: the point is currently highlighted
    background: the point is background information
*/

// for JS
export const flagSelected = 1;
export const flagBackground = 2;
export const flagHighlight = 4;

// for GLSL
export const glPointFlags = `

  const float flagSelected = 1.;
  const float flagBackground = 2.;
  const float flagHighlight = 4.;

  bool isLowBitSet(float f) {
    f = mod(f, 2.);
    return (f > 0.9 && f <= 1.1);
  }

  float shiftRightOne(float f) {
    return floor(f / 2.);
  }

  void getFlags(in float flag,
                out bool isBackground,
                out bool isSelected,
                out bool isHighlight) {
    isSelected = isLowBitSet(flag);
    flag = shiftRightOne(flag);
    isBackground = isLowBitSet(flag);
    flag = shiftRightOne(flag);
    isHighlight = isLowBitSet(flag);
  }

`;

/*
Point Size:
  Calculate point size for scatter plot based upon pseudo density.

  Linear scaling of point size, clamped to [minSize, maxSize],
  between two anchor points based on empirical testing.

  Anchor points:
    - 5M points on a 500x500 canvas: density=20 -> 1.0 (minimum visible)
    - 1000 points on a 1440x1440 canvas: density~=0.0005 -> 5.0

  The domain is pseudo density (numPoints / minViewportDimension^2).
  The range is WebGL point size.

  The minimum is 1.0 because WebGL implementations clamp gl_PointSize
  to ALIASED_POINT_SIZE_RANGE (typically [1, max]). Sub-pixel point
  sizes cause rendering inconsistencies across GPU drivers and make
  selection state invisible. See issue #2709.
*/

// configuration - extended domain to support datasets up to ~10M cells
const domain = [5000000 / (500 * 500), 1000 / (1440 * 1440)];
const range = [1.0, 5];

// derived from configuration
const scale = (range[1] - range[0]) / (domain[1] - domain[0]);
const offset = scale * -domain[0] + range[0];

export const glPointSize = `
  float pointSize(float nPoints, float minViewportDimension, bool isSelected, bool isHighlight) {
    float density = nPoints / (minViewportDimension * minViewportDimension);
    float pointSize = (${scale.toFixed(4)}*density) + ${offset.toFixed(4)};
    pointSize = clamp(pointSize, 
      ${range[0].toFixed(4)}, 
      ${range[1].toFixed(4)});

    if (isHighlight) return max(2. * pointSize, 2.0);
    if (isSelected) return max(pointSize, 1.0);
    return max(pointSize / 3., 0.5);
  }
`;

/*
Point Alpha:
  Provide alpha-based differentiation between selected and unselected
  cells, so the selection state remains visible even when point sizes
  are at their minimum. This ensures large datasets (4M+ cells)
  are visually distinguishable. See issue #2709.
*/
export const glPointAlpha = `
  float pointAlpha(bool isBackground, bool isSelected, bool isHighlight) {
    if (isHighlight) return 1.0;
    if (isBackground) return 0.9;
    if (!isSelected) return 0.3;
    return 1.0;
  }
`;
