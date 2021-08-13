/*
Private helper functions related to schema
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function _getColumnSchema(schema: any, field: any, col: any) {
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
      throw new Error(`unknown field name: ${field}`);
  }
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function _isIndex(schema: any, field: any, col: any): boolean {
  const index = schema.annotations?.[field].index;
  return index && index === col;
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export function _getColumnDimensionNames(schema: any, field: any, col: any) {
  /*
		field/col may be an alias for multiple columns. Currently used to map ND 
		values to 1D dataframe columns for embeddings/layout. Signified by the presence
		of the "dims" value in the schema.
		*/
  const colSchema = _getColumnSchema(schema, field, col);
  if (!colSchema) {
    return undefined;
  }
  return colSchema.dims || [col];
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'schema' implicitly has an 'any' type.
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function _schemaColumns(schema, field) {
  switch (field) {
    case "obs":
      return Object.keys(schema.annotations.obsByName);
    case "var":
      return Object.keys(schema.annotations.varByName);
    case "emb":
      return Object.keys(schema.layout.obsByName);
    default:
      throw new Error(`unknown field name: ${field}`);
  }
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'schema' implicitly has an 'any' type.
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function _getWritableColumns(schema, field) {
  if (field !== "obs") return [];
  return (
    schema.annotations.obs.columns
      // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'v' implicitly has an 'any' type.
      .filter((v) => v.writable)
      // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'v' implicitly has an 'any' type.
      .map((v) => v.name)
  );
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'schema' implicitly has an 'any' type.
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function _isContinuousType(schema) {
  const { type } = schema;
  return !(type === "string" || type === "boolean" || type === "categorical");
}
