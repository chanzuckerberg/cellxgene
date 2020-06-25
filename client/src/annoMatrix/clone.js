export default function clone(orig) {
  return Object.assign(Object.create(Object.getPrototypeOf(orig)), orig);
}
