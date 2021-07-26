/*
Various type and schema related helper functions.
*/

/*
Utility function to test for a typed array
*/
export function isTypedArray(x: any) {
  return (
    ArrayBuffer.isView(x) &&
    Object.prototype.toString.call(x) !== "[object DataView]"
  );
}

/*
Test for float typed array, ie, Float32TypedArray or Float64TypedArray
*/
export function isFpTypedArray(x: any) {
  let constructor;
  const isFloatArray =
    x &&
    ({ constructor } = x) &&
    (constructor === Float32Array || constructor === Float64Array);
  return isFloatArray;
}

export function isArrayOrTypedArray(x: any) {
  return Array.isArray(x) || isTypedArray(x);
}
