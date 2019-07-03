/*
Helpers for schema management
*/
import _ from "lodash";

import fromEntries from "../fromEntries";

/*
System wide schema assumptions:
	- schema and data wil be consistent (eg, for user-created annotations)
	- schema will be internally self-consistent (eg, index matches columns)
	- world & universe schema are same - only data is subset
*/

export function indexEntireSchema(schema) {
	/* Index schema for ease of use */
	schema.annotations.obsByName = fromEntries(
		schema.annotations.obs.columns.map(v => [v.name, v])
	);
	schema.annotations.varByName = fromEntries(
		schema.annotations.var.columns.map(v => [v.name, v])
	);
	schema.layout.obsByName = fromEntries(
		schema.layout.obs.map(v => [v.name, v])
	);
	schema.layout.varByName = fromEntries(
		schema.layout.var.map(v => [v.name, v])
	);

	return schema;
}

function _copy(schema) {
	/* redux copy conventions - WARNING, only for modifyign obs annotations */
	return {
		...schema,
		annotations: {
			...schema.annotations,
			obs: _.cloneDeep(schema.annotations.obs)
		}
	};
}

function _reindex(schema) {
	/* reindex obs annotations ONLY */
	schema.annotations.obsByName = fromEntries(
		schema.annotations.obs.columns.map(v => [v.name, v])
	);
	return schema;
}

export function removeObsAnnoColumn(schema, name) {
	const newSchema = _copy(schema);
	newSchema.annotations.obs.columns = schema.annotations.obs.columns.filter(
		v => v.name !== name
	);
	return _reindex(newSchema);
}

export function addObsAnnoColumn(schema, name, defn) {
	const newSchema = _copy(schema);
	newSchema.annotations.obs.columns.push(defn);
	return _reindex(newSchema);
}

export function removeObsAnnoCategory(schema, name, category) {
	/* remove a category from a categorical annotation */
	const categories = schema.annotations.obsByName[name]?.categories;
	if (!categories)
		throw new Error("column does not exist or is not categorical");

	const idx = categories.indexOf(category);
	if (idx === -1) throw new Error("category does not exist");

	const newSchema = _reindex(_copy(schema));

	/* remove category */
	newSchema.annotations.obsByName[name].categories.splice(idx, 1);
	return newSchema;
}

export function addObsAnnoCategory(schema, name, category) {
	/* add a category to a categorical annotation */
	const categories = schema.annotations.obsByName[name]?.categories;
	if (!categories)
		throw new Error("column does not exist or is not categorical");

	const idx = categories.indexOf(category);
	if (idx !== -1) throw new Error("category already exists");

	const newSchema = _reindex(_copy(schema));

	/* remove category */
	newSchema.annotations.obsByName[name].categories.push(category);
	return newSchema;
}
