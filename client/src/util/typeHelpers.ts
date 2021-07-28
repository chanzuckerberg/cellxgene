/*
Various type and schema related helper functions.
*/

/*
Utility function to test for a typed array
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function isTypedArray(x: any) {
  return (
    ArrayBuffer.isView(x) &&
    Object.prototype.toString.call(x) !== "[object DataView]"
  );
}

/*
Test for float typed array, ie, Float32TypedArray or Float64TypedArray
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function isFpTypedArray(x: any) {
  let constructor;
  const isFloatArray =
    x &&
    ({ constructor } = x) &&
    (constructor === Float32Array || constructor === Float64Array);
  return isFloatArray;
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function isArrayOrTypedArray(x: any) {
  return Array.isArray(x) || isTypedArray(x);
}
