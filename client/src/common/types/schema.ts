export type Category = number | string | boolean;

export interface AnnotationColumnSchema {
  categories?: Category[];
  name: string;
  type: "string" | "float32" | "int32" | "categorical" | "boolean";
  writable: boolean;
}

export interface XMatrixSchema {
  nObs: number;
  nVar: number;
  // TODO(thuang): Not sure what other types are available
  type: "float32";
}

export interface EmbeddingSchema {
  dims: string[];
  name: string;
  // TODO(thuang): Not sure what other types are available
  type: "float32";
}
interface RawLayoutSchema {
  obs: EmbeddingSchema[];
  var?: EmbeddingSchema[];
}

interface RawAnnotationsSchema {
  obs: {
    columns: AnnotationColumnSchema[];
    index: string;
  };
  var: {
    columns: AnnotationColumnSchema[];
    index: string;
  };
}

export interface RawSchema {
  annotations: RawAnnotationsSchema;
  dataframe: XMatrixSchema;
  layout: RawLayoutSchema;
}

interface AnnotationsSchema extends RawAnnotationsSchema {
  obsByName: { [name: string]: AnnotationColumnSchema };
  varByName: { [name: string]: AnnotationColumnSchema };
}

interface LayoutSchema extends RawLayoutSchema {
  obsByName: { [name: string]: EmbeddingSchema };
  varByName: { [name: string]: EmbeddingSchema };
}

export interface Schema extends RawSchema {
  annotations: AnnotationsSchema;
  layout: LayoutSchema;
}

/**
 * Sub-schema objects describing the schema for a primitive Array or Matrix in one of the fields.
 */
export type ArraySchema =
  | AnnotationColumnSchema
  | EmbeddingSchema
  | XMatrixSchema;

/**
 * Set of data / metadata objects that must be specified in a CXG.
 */
export enum Field {
  "obs" = "obs",
  "var" = "var",
  "emb" = "emb",
  "X" = "X",
}
