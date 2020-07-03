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

function compare(a, b) {
  const diff = a.priority - b.priority;
  if (diff) return diff;
  return a.order - b.order;
}

export default class PromiseLimit {
  constructor(maxConcurrency = 5) {
    this.queue = new TinyQueue([], compare);
    this.maxConcurrency = maxConcurrency;
    this.pending = 0;
    this.insertCounter = 0;
  }

  priorityAdd(p, fn, ...args) {
    // p - numermic priority (lower first)
    // fn - must return a promise
    // args - will be passed to fn
    return this._push(p, fn, args);
  }

  add(fn, ...args) {
    // fn - must return a promise
    // args - will be passed to fn
    return this._push(0, fn, args);
  }

  /**
  Private below
  **/

  _push(priority, fn, args) {
    const order = this.insertCount;
    this.insertCount += 1;
    return new Promise((resolve, reject) => {
      this.queue.push({ priority, order, fn, args, resolve, reject });
      this._resolveNext(false);
    });
  }

  _resolveNext = (completed = true) => {
    if (completed) this.pending -= 1;

    while (this.queue.length > 0 && this.pending < this.maxConcurrency) {
      const task = this.queue.pop(); // order of insertion
      this.pending += 1;
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
