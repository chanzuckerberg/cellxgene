/*
Shallow clone an object, correctly handling prototype
*/
export default function _shallowClone(orig) {
  return Object.assign(Object.create(Object.getPrototypeOf(orig)), orig);
}
