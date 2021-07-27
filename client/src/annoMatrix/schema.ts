/*
Private helper functions related to schema
*/
import catLabelSort from "../util/catLabelSort";
import { unassignedCategoryLabel } from "../globals";

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
export function _isContinuousType(schema) {
  const { type } = schema;
  return !(type === "string" || type === "boolean" || type === "categorical");
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'colSchema' implicitly has an 'any' type... Remove this comment to see the full error message
export function _normalizeCategoricalSchema(colSchema, col) {
  const { type, writable } = colSchema;
  if (
    type === "string" ||
    type === "boolean" ||
    type === "categorical" ||
    writable
  ) {
    const categorySet = new Set(
      col.summarizeCategorical().categories.concat(colSchema.categories ?? [])
    );
    if (writable && !categorySet.has(unassignedCategoryLabel)) {
      categorySet.add(unassignedCategoryLabel);
    }
    colSchema.categories = Array.from(categorySet);
  }

  if (colSchema.categories) {
    colSchema.categories = catLabelSort(writable, colSchema.categories);
  }
  return colSchema;
}
