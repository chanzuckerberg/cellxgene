/*
Helper functions for user-editable nnotations state management.
See also reducers/annotations.js
*/
import { unassignedCategoryLabel } from "../../globals";
import * as SchemaHelpers from "./schemaHelpers";

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
	return SchemaHelpers.removeObsAnnoColumn(schema, name);
}

export function addObsAnnoSchema(schema, name, colSchema) {
	/*
	add a categorical type to the obs annotation schema 
	*/

	/* collision detection */
	if (schema.annotations.obs.columns.some(v => v.name === name))
		throw Error("annotations may not contain duplicate category names");
	if (name !== colSchema.name) throw Error("column schema does not match");
	return SchemaHelpers.addObsAnnoColumn(schema, name, colSchema);
}

export function removeObsAnnoCategory(schema, name, category) {
	/* don't allow deletion of unassigned category on writable annotations */

	if (!_isUserAnnotation(schema, name))
		throw new Error("unable to modify read-only schema");
	if (category === unassignedCategoryLabel)
		throw new Error("may not remove unassigned category lable");

	return SchemaHelpers.removeObsAnnoCategory(schema, name, category);
}

export function addObsAnnoCategory(schema, name, category) {
	if (!_isUserAnnotation(schema, name))
		throw new Error("unable to modify read-only schema");

	return SchemaHelpers.addObsAnnoCategory(schema, name, category);
}

export function setLabelByValue(df, colName, fromLabel, toLabel) {
	/* 
	in the dataframe column `colName`, set any value of `fromLabel` to `toLabel`
	*/
	const keys = df.colIndex.keys();
	const ndf = df.mapColumns((col, colIdx) => {
		if (colName !== keys[colIdx]) return col;

		/* clone data and return it. */
		const newCol = col.slice();
		for (let i = 0, l = newCol.length; i < l; i += 1) {
			if (newCol[i] === fromLabel) newCol[i] = toLabel;
		}
		return newCol;
	});
	return ndf;
}

export function setLabelByMask(df, colName, mask, label) {
	/* 
	in the dataframe column `colName`, set the masked rows to 'label'
	*/
	const keys = df.colIndex.keys();
	const ndf = df.mapColumns((col, colIdx) => {
		if (colName !== keys[colIdx]) return col;

		/* clone data and return it. */
		const newCol = col.slice();
		for (let i = 0, l = newCol.length; i < l; i += 1) {
			if (mask[i]) newCol[i] = label;
		}
		return newCol;
	});
	return ndf;
}

export function worldToUniverseMask(worldMask, worldObsAnnotations, nObs) {
	const mask = new Uint8Array(nObs);
	const { rowIndex } = worldObsAnnotations;

	for (let i = 0, l = worldMask.length; i < l; i += 1) {
		if (worldMask[i]) {
			const label = rowIndex.getLabel(i);
			mask[label] = 1;
		}
	}
	return mask;
}
