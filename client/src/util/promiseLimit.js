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
*/
export default class PromiseLimit {
  constructor(maxConcurrency) {
    this.queue = new Set();
    this.maxConcurrency = maxConcurrency;
    this.pending = 0;
  }

  add(fn, ...args) {
    // fn - must return a promise
    // args - will be passed to fn
    return new Promise((resolve, reject) => {
      this.queue.add({ fn, args, resolve, reject });
      this._resolveNext(false);
    });
  }

  _resolveNext = (completed = true) => {
    if (completed) this.pending -= 1;

    while (this.queue.size > 0 && this.pending < this.maxConcurrency) {
      const task = this.queue.values().next().value; // order of insertion
      this.pending += 1;
      this.queue.delete(task);
      const { resolve, reject, fn, args } = task;

      try {
        const result = fn(...args);
        result.then(this._resolveNext, this._resolveNext);
        result.then(resolve, reject);
      } catch (err) {
        reject(err);
      }
    }
  };
}
