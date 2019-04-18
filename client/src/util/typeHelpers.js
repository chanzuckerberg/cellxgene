/*
Various type and schema related helper functions.
*/

/*
Utility function to test for a typed array
*/
export function isTypedArray(x) {
	return (
		ArrayBuffer.isView(x) &&
		Object.prototype.toString.call(x) !== "[object DataView]"
	);
}

/*
Test for float typed array, ie, Float32TypedArray or Float64TypedArray
*/
export function isFpTypedArray(x) {
	let ctor;
	const ret =
		x &&
		(ctor = x.constructor) &&
		(ctor === Float32Array || ctor === Float64Array);
	return ret;
}

export function isArrayOrTypedArray(x) {
	return Array.isArray(x) || isTypedArray(x);
}
