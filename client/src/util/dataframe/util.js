/*
Private utility code for dataframe
*/

export { isTypedArray, isArrayOrTypedArray } from "../typeHelpers";

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
