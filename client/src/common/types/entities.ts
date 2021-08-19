// If a globally shared type or interface doesn't have a clear owner, put it here

export type Category = number | string | boolean;

export interface AnnotationColumn {
  categories?: Category[];
  name: string;
  type: "string" | "float32" | "int32" | "categorical" | "boolean";
  writable: boolean;
}

interface DataFrame {
  nObs: number;
  nVar: number;
  // TODO(thuang): Not sure what other types are available
  type: "float32";
}

export interface LayoutColumn {
  dims: string[];
  name: string;
  // TODO(thuang): Not sure what other types are available
  type: "float32";
}
interface RawLayout {
  obs: LayoutColumn[];
  var?: LayoutColumn[];
}

interface RawAnnotations {
  obs: {
    columns: AnnotationColumn[];
    index: string;
  };
  var: {
    columns: AnnotationColumn[];
    index: string;
  };
}

export interface RawSchema {
  annotations: RawAnnotations;
  dataframe: DataFrame;
  layout: RawLayout;
}

interface Annotations extends RawAnnotations {
  obsByName: { [name: string]: AnnotationColumn };
  varByName: { [name: string]: AnnotationColumn };
}

interface Layout extends RawLayout {
  obsByName: { [name: string]: LayoutColumn };
  varByName: { [name: string]: LayoutColumn };
}

export interface Schema extends RawSchema {
  annotations: Annotations;
  layout: Layout;
}
