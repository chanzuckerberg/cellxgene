/*
AnnoMatrix -- Annotated Matrix exported interface

Public API is defined in:

	annoMatrix.js
	viewCreators.js
	crossfilter.js

*/

export { default as AnnoMatrixLoader } from "./loader";
export * from "./viewCreators";
export { default as AnnoMatrixObsCrossfilter } from "./crossfilter";
export { default as gcMiddleware } from "./middleware";
