/*
Utility code for our WebGL shaders
*/

/*
Point flags - used in graph & scatter plots

We want a bitmask-like flag structure, but due to webgl limitations
must emulate it with floats.
*/

// for JS
export const flagSelected = 1;
export const flagNaN = 2;
export const flagHighlight = 4;

// for GLSL
export const glPointFlags = `

  const float flagSelected = 1.;
  const float flagNaN = 2.;
  const float flagHighlight = 4.;

  bool isLowBitSet(float f) {
    f = mod(f, 2.);
    return (f > 0.9 && f <= 1.1);
  }

  float shiftRightOne(float f) {
    return floor(f / 2.);
  }

  void getFlags(in float flag,
                out bool isNaN,
                out bool isSelected,
                out bool isHighlight) {
    isSelected = isLowBitSet(flag);
    flag = shiftRightOne(flag);
    isNaN = isLowBitSet(flag);
    flag = shiftRightOne(flag);
    isHighlight = isLowBitSet(flag);
  }

`;
