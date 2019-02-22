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

export function callOnceLazy(f) {
  let value;
  let calledOnce = false;
  const result = function result(...args) {
    if (!calledOnce) {
      value = f(...args);
      calledOnce = true;
    }
    return value;
  };

  return result;
}
