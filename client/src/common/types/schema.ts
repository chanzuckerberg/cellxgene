type Category = number | string | boolean;

export interface AnnotationColumnSchema {
  categories?: Category[];
  name: string;
  type: "string" | "float32" | "int32" | "categorical" | "boolean";
  writable: boolean;
}

interface DataFrameSchema {
  nObs: number;
  nVar: number;
  // TODO(thuang): Not sure what other types are available
  type: "float32";
}

export interface LayoutColumnSchema {
  dims: string[];
  name: string;
  // TODO(thuang): Not sure what other types are available
  type: "float32";
}
interface RawLayoutSchema {
  obs: LayoutColumnSchema[];
  var?: LayoutColumnSchema[];
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
  dataframe: DataFrameSchema;
  layout: RawLayoutSchema;
}

interface AnnotationsSchema extends RawAnnotationsSchema {
  obsByName: { [name: string]: AnnotationColumnSchema };
  varByName: { [name: string]: AnnotationColumnSchema };
}

interface LayoutSchema extends RawLayoutSchema {
  obsByName: { [name: string]: LayoutColumnSchema };
  varByName: { [name: string]: LayoutColumnSchema };
}

export interface Schema extends RawSchema {
  annotations: AnnotationsSchema;
  layout: LayoutSchema;
}
