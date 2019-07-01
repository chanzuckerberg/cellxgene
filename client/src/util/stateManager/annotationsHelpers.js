/*
Helper functions for user-editable nnotations state management.
See also reducers/annotations.js
*/

import { removeObsAnnoColumn, addObsAnnoColumn } from "./schemaHelpers";

/*
There are a number of state constraints assumed throughout the
application:
	- all obs annotations are in {world|universe}.obsAnnotations,
		regardless of whether or not they are user editable.
	- the {world|universe}.schema is always up to date and matches
		the data
	- the schema flag `writable` correctly indicates whether
	  the annotation is editable/mutable.
*/

function _isUserAnnotation(schema, name) {
	return schema.annotations.obsByName[name]?.writable;
}

export function isUserAnnotation(worldOrUniverse, name) {
	return _isUserAnnotation(worldOrUniverse.schema, name);
}

export function removeObsAnnoSchema(schema, name) {
	/*
	remove named annotation from obs annotation schema
	*/

	/* only remove if it exists and is a user annotation */
	if (!_isUserAnnotation(schema, name))
		throw new Error("removing non-user-defined schema");
	return removeObsAnnoColumn(schema, name);
}

export function addObsAnnoSchema(schema, name, colSchema) {
	/*
	add a categorical type to the obs annotation schema 
	*/

	/* collision detection */
	if (schema.annotations.obs.columns.some(v => v.name === name))
		throw Error("annotations may not contain duplicate category names");
	if (name !== colSchema.name) throw Error("column schema does not match");
	return addObsAnnoColumn(schema, name, colSchema);
}
