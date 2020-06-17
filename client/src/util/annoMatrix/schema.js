/*
Helper functions related to schema
*/

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
