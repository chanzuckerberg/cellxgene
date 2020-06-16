/*
Helper functions related to schema
*/

import fromEntries from "../fromEntries";

export function getColumnSchema(schema, field, col) {
  /* look up the column definition */
  switch (field) {
    case "obs":
      if (typeof col === "object")
        throw new Error("unable to get column schema by query");
      return schema.annotations.obsByName[col];
    case "var":
      if (typeof col === "object")
        throw new Error("unable to get column schema by query");
      return schema.annotations.varByName[col];
    case "emb":
      if (typeof col === "object")
        throw new Error("unable to get column schema by query");
      return schema.layout.obsByName[col];
    case "X":
      return schema.dataframe;
    default:
      throw new Error("unknown field name");
  }
}

export function getColumnDimensionNames(schema, field, col) {
  /*
		field/col may be an alias for multiple columns. Currently used to map ND 
		values to 1D dataframe columns for embeddings/layout. Signfied by the presence
		of the "dims" value in the schema.
		*/
  const colSchema = getColumnSchema(schema, field, col);
  if (!colSchema) {
    throw new Error("unknown column name");
  }
  return colSchema.dims || [col];
}

export function schemaColumns(schema, field) {
  switch (field) {
    case "obs":
      return Object.keys(schema.annotations.obsByName);
    case "var":
      return Object.keys(schema.annotations.varByName);
    case "emb":
      return Object.keys(schema.layout.obsByName);
    default:
      throw new Error("unknown field name");
  }
}

export function isContinuousType(schema) {
  const { type } = schema;
  return !(type === "string" || type === "boolean" || type === "categorical");
}

export function indexEntireSchema(schema) {
  /* Index schema for ease of use */
  schema.annotations.obsByName = fromEntries(
    schema.annotations?.obs?.columns?.map((v) => [v.name, v]) ?? []
  );
  schema.annotations.varByName = fromEntries(
    schema.annotations?.var?.columns?.map((v) => [v.name, v]) ?? []
  );
  schema.layout.obsByName = fromEntries(
    schema.layout?.obs?.map((v) => [v.name, v]) ?? []
  );
  schema.layout.varByName = fromEntries(
    schema.layout?.var?.map((v) => [v.name, v]) ?? []
  );

  return schema;
}

// function fromEntries(arr) {
/*
	Similar to Object.fromEntries, but only handles array.
	This could be replaced with the standard fucnction once it
	is widely available.   As of 3/20/2019, it has not yet
	been released in the Chrome stable channel.
	*/
// 	const obj = {};
// 	for (let i = 0, l = arr.length; i < l; i += 1) {
// 		obj[arr[i][0]] = arr[i][1];
// 	}
// 	return obj;
// }
