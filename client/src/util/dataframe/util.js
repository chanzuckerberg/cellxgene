/*
Private utility code for dataframe
*/

export function isTypedArray(x) {
  return (
    ArrayBuffer.isView(x) &&
    Object.prototype.toString.call(x) !== "[object DataView]"
  );
}

export function isArrayOrTypedArray(x) {
  return Array.isArray(x) || isTypedArray(x);
}
