export default function _clone(orig) {
  return Object.assign(Object.create(Object.getPrototypeOf(orig)), orig);
}
