import { NumericArray } from "../common/types/entities";

/*
clip - clip all values in a Array or TypedArray, IN PLACE.

Values in array are clipped if less than `lower` or greater than `upper`.

If `setTo` is undefined, values less than `lower` will be set to `lower`,
and values greater than `upper` will be set to `upper`.

If `setTo` is not undefined, values outside the [lower, upper] range will be set to
`setTo`.

*/
export default function clip(
  arr: NumericArray,
  lower: number,
  upper: number,
  setTo?: number
): NumericArray {
  const lowerSet = setTo === undefined ? lower : setTo;
  const upperSet = setTo === undefined ? upper : setTo;
  for (let i = 0, l = arr.length; i < l; i += 1) {
    const v = arr[i];
    if (v < lower) {
      arr[i] = lowerSet;
    } else if (v > upper) {
      arr[i] = upperSet;
    }
  }
  return arr;
}
