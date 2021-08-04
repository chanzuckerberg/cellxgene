/*
Array range creation

range(start, stop, step) -> Array
	This is identical to https://docs.python.org/3/library/functions.html#func-range
	Returns new array filled with a range of numbers.

	Usage:

		range(stop) - start defaults to zero, step defaults to 1
		range(start, stop, [step]) - step defaults to 1

	Examples:
		range(3) -> [0, 1, 2]
		range(1, 3) -> [1, 2]
		range(1, 5, 2) -> [1, 3]


rangeFill(array, start, step) -> array
	Fill entire array with values, from start, by step.  Returns first array.
	start defaults to zero, step defaults to one.

*/

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'arr' implicitly has an 'any' type.
function _doFill(arr, start, step, count) {
  for (let idx = 0, val = start; idx < count; idx += 1, val += step) {
    arr[idx] = val;
  }
  return arr;
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'arr' implicitly has an 'any' type.
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function rangeFill(arr, start = 0, step = 1) {
  return _doFill(arr, start, step, arr.length);
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'start' implicitly has an 'any' type.
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function range(start, stop, step) {
  if (start === undefined) return [];
  if (stop === undefined) {
    stop = start;
    start = 0;
  }
  step = step || 1; // catch undefind and zero
  const len = Math.max(Math.ceil((stop - start) / step), 0);
  return _doFill(new Array(len), start, step, len);
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'start' implicitly has an 'any' type.
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function linspace(start, stop, nsteps) {
  // @ts-expect-error ts-migrate(2363) FIXME: The right-hand side of an arithmetic operation mus... Remove this comment to see the full error message
  const delta = (stop - start) / (nsteps - 1).toFixed();
  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'i' implicitly has an 'any' type.
  return range(0, nsteps, 1).map((i) => start + i * delta);
}
