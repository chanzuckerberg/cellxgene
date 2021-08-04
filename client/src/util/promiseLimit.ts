/*
This class provides a means to resolve Promises asynchronously,
and limit the number concurrently executed.  For example, if
you want to involve a large number of API endpoints using fetch(),
but want to limit the number simultaneously outstanding.

Example usage:

const plimit = new PromiseLimit(2);
return Promise.all([
	plimit.add(() => fetch('/foo')),
	plimit.add(() => fetch('/bar')),
	plimit.add(() => fetch('/baz'))
])

if you want a priority queue based implementation, just
use priorityAdd() instead of add():

const plimit = new PromiseLimit(2);
return Promise.all([
  plimit.priorityAdd(0, () => fetch('/foo')),
  plimit.priorityAdd(10, () => fetch('/bar')),
  plimit.priorityAdd(-1, () => fetch('/baz'))
])

Priority is a numeric value.  Lower first.  Stable ordering.
*/

import TinyQueue from "tinyqueue";

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'a' implicitly has an 'any' type.
function compare(a, b) {
  const diff = a.priority - b.priority;
  if (diff) return diff;
  return a.order - b.order;
}

export default class PromiseLimit {
  constructor(maxConcurrency = 5) {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (this as any).queue = new TinyQueue([], compare);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (this as any).maxConcurrency = maxConcurrency;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (this as any).pending = 0;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (this as any).insertCounter = 0;
  }

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'p' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  priorityAdd(p, fn, ...args) {
    // p - numermic priority (lower first)
    // fn - must return a promise
    // args - will be passed to fn
    return this._push(p, fn, args);
  }

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'fn' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  add(fn, ...args) {
    // fn - must return a promise
    // args - will be passed to fn
    return this._push(0, fn, args);
  }

  /**
  Private below
  **/

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'priority' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  _push(priority, fn, args) {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const order = (this as any).insertCount;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (this as any).insertCount += 1;
    return new Promise((resolve, reject) => {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (this as any).queue.push({ priority, order, fn, args, resolve, reject });
      this._resolveNext(false);
    });
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  _resolveNext = (completed = true) => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    if (completed) (this as any).pending -= 1;

    while (
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (this as any).queue.length > 0 &&
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (this as any).pending < (this as any).maxConcurrency
    ) {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      const task = (this as any).queue.pop(); // order of insertion
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (this as any).pending += 1;
      const { resolve, reject, fn, args } = task;

      try {
        const result = fn(...args);
        result.then(
          () => this._resolveNext(true),
          () => this._resolveNext(true)
        );
        result.then(resolve, reject);
      } catch (err) {
        reject(err);
      }
    }
  };
}
