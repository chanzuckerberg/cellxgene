/*
Private helper functions related to schema
*/
import {
  AnnotationColumnSchema,
  ArraySchema,
  Field,
  Schema,
} from "../common/types/schema";
import { LabelArray, LabelType } from "../util/dataframe/types";

export function _getColumnSchema(
  schema: Schema,
  field: Field,
  col: LabelType
): ArraySchema {
  /* look up the column definition */
  switch (field) {
    case Field.obs:
      if (typeof col === "object")
        throw new Error("unable to get column schema by query");
      return schema.annotations.obsByName[col];
    case Field.var:
      if (typeof col === "object")
        throw new Error("unable to get column schema by query");
      return schema.annotations.varByName[col];
    case Field.emb:
      if (typeof col === "object")
        throw new Error("unable to get column schema by query");
      return schema.layout.obsByName[col];
    case Field.X:
      return schema.dataframe;
    default:
      throw new Error(`unknown field name: ${field}`);
  }
}

export function _isIndex(
  schema: Schema,
  field: Field.obs | Field.var,
  col: LabelType
): boolean {
  const index = schema.annotations?.[field].index;
  return !!(index && index === col);
}

export function _getColumnDimensionNames(
  schema: Schema,
  field: Field,
  col: LabelType
): LabelArray | undefined {
  /*
		field/col may be an alias for multiple columns. Currently used to map ND 
		values to 1D dataframe columns for embeddings/layout. Signified by the presence
		of the "dims" value in the schema.
		*/
  const colSchema = _getColumnSchema(schema, field, col);
  if (!colSchema) {
    return undefined;
  }
  if ("dims" in colSchema) {
    return colSchema.dims;
  }
  return [col];
}

export function _schemaColumns(schema: Schema, field: Field): string[] {
  switch (field) {
    case Field.obs:
      return Object.keys(schema.annotations.obsByName);
    case Field.var:
      return Object.keys(schema.annotations.varByName);
    case Field.emb:
      return Object.keys(schema.layout.obsByName);
    default:
      throw new Error(`unknown field name: ${field}`);
  }
}

export function _getWritableColumns(schema: Schema, field: Field): string[] {
  if (field !== Field.obs) return [];
  return schema.annotations.obs.columns
    .filter((v: AnnotationColumnSchema) => v.writable)
    .map((v: AnnotationColumnSchema) => v.name);
}

export function _isContinuousType(schema: ArraySchema): boolean {
  const { type } = schema;
  return !(type === "string" || type === "boolean" || type === "categorical");
}
