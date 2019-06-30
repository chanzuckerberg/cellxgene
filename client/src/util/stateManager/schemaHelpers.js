/*
Helpers for schema management
*/
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

export function removeObsAnnoColumn(schema, name) {
	/* redux copy conventions */
	const newSchema = {
		...schema,
		annotations: {
			...schema.annotations,
			obs: {
				...schema.annotations.obs,
				columns: schema.annotations.obs.columns.filter(v => v.name !== name)
			}
		}
	};

	/* reindex */
	newSchema.annotations.obsByName = fromEntries(
		newSchema.annotations.obs.columns.map(v => [v.name, v])
	);

	return newSchema;
}

export function addObsAnnoColumn(schema, name, defn) {
	const newSchema = {
		...schema,
		annotations: {
			...schema.annotations,
			obs: {
				...schema.annotations.obs,
				columns: schema.annotations.obs.columns.concat([defn])
			}
		}
	};

	/* reindex */
	newSchema.annotations.obsByName = fromEntries(
		newSchema.annotations.obs.columns.map(v => [v.name, v])
	);

	return newSchema;
}
