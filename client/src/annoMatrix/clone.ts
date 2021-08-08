/*
Shallow clone an object, correctly handling prototype
*/

export default function _shallowClone<T>(orig: T): T {
  return Object.assign(Object.create(Object.getPrototypeOf(orig)), orig);
}
