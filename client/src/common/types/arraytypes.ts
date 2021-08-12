/**
 * Utility type and interface definitions.
 */

/**
 * TypedArrays that can be assigned to a number.
 */
export type TypedArray =
  | Int8Array
  | Uint8Array
  | Int16Array
  | Uint16Array
  | Int32Array
  | Uint32Array
  | Float32Array
  | Float64Array;

export type UnsignedTypedArray = Uint8Array | Uint16Array | Uint32Array;
export type FloatTypedArray = Float32Array | Float64Array;

export type TypedArrayConstructor =
  | Int8ArrayConstructor
  | Uint8ArrayConstructor
  | Int16ArrayConstructor
  | Uint16ArrayConstructor
  | Int32ArrayConstructor
  | Uint32ArrayConstructor
  | Float32ArrayConstructor
  | Float64ArrayConstructor;

export type AnyArray = Array<unknown> | TypedArray;

export type NumberArray = Array<number> | TypedArray;

export type Int8 = Int8Array[0];
export type Uint8 = Uint8Array[0];
export type Int16 = Int16Array[0];
export type Uint16 = Uint16Array[0];
export type Int32 = Int32Array[0];
export type Uint32 = Uint32Array[0];
export type Float32 = Float32Array[0];
export type Float64 = Float64Array[0];

/**
 * Test if the parameter is a TypedArray.
 * @param tbd - value to be tested
 * @returns true if `tbd` is a TypedArray, false if not.
 */
export function isTypedArray(tbd: unknown): tbd is TypedArray {
  return (
    ArrayBuffer.isView(tbd) &&
    Object.prototype.toString.call(tbd) !== "[object DataView]"
  );
}

/**
 * Test if the paramter is a float TypedArray
 * @param tbd - value to be tested
 * @returns - true if `tbd` is a float typed array.
 */
export function isFloatTypedArray(tbd: unknown): tbd is FloatTypedArray {
  return tbd instanceof Float32Array || tbd instanceof Float64Array;
}

/**
 * Test if the paramter is a float TypedArray
 * @param tbd - value to be tested
 * @returns - true if `tbd` is a float typed array.
 */
export function isUnsignedTypedArray(tbd: unknown): tbd is UnsignedTypedArray {
  return (
    tbd instanceof Uint8Array ||
    tbd instanceof Uint16Array ||
    tbd instanceof Uint32Array
  );
}

/**
 * Test if the parameter is a TypedArray or Array
 * @param tbd - value to be tested
 * @returns - true if `tbd` is a TypedArray or Array
 */
export function isAnyArray(tbd: unknown): tbd is AnyArray {
  return Array.isArray(tbd) || isTypedArray(tbd);
}
