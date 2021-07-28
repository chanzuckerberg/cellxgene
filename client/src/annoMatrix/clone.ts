/*
Shallow clone an object, correctly handling prototype
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export default function _shallowClone(orig: any) {
  return Object.assign(Object.create(Object.getPrototypeOf(orig)), orig);
}
