import { TypedArray } from "../../common/types/arraytypes";

export type LabelType = number | string;

type CommonProps<A, B> = {
  [K in keyof A & keyof B]: A[K] | B[K];
};
export type GenericLabelArray<T> = CommonProps<Array<T>, Int32Array>;

export type RowLabelArray = GenericLabelArray<number>;
export type ColLabelArray = GenericLabelArray<number | string>;
export type LabelArray = GenericLabelArray<number | string>;

export type OffsetType = number;
export type OffsetArray =
  | Int8Array
  | Uint8Array
  | Int16Array
  | Uint16Array
  | Int32Array
  | Uint32Array
  | number[];

export type ContinuousColumnSummary = {
  categorical: false;
  min: number;
  max: number;
  nan: number;
  pinf: number;
  ninf: number;
  percentiles: number[];
};

export type CategoricalColumnSummary = {
  categorical: true;
  categories: (number | string | boolean)[];
  categoryCounts: Map<number | string | boolean, number>;
  numCategories: number;
};

export type ColumnSummary = ContinuousColumnSummary | CategoricalColumnSummary;

export type ContinuousHistogram = number[];
export type ContinuousHistogramBy = Map<DataframeValue, ContinuousHistogram>;
export type CategoricalHistogram = Map<DataframeValue, number>;
export type CategoricalHistogramBy = Map<DataframeValue, CategoricalHistogram>;

export type DataframeValue = number | string | boolean;

export type DataframeColumnarArray = DataframeValue[] | TypedArray;

export type DataframeColumnGetter = (
  label: LabelType
) => DataframeValue | undefined;

export interface DataframeColumn extends DataframeColumnGetter {
  readonly __id: string;
  isContinuous: boolean;

  asArray: () => DataframeColumnarArray;

  summarizeContinuous: () => ContinuousColumnSummary;
  summarizeCategorical: () => CategoricalColumnSummary;

  histogramContinuous: (
    bins: number,
    domain: [number, number]
  ) => ContinuousHistogram;
  histogramContinuousBy: (
    bins: number,
    domain: [number, number],
    by: DataframeColumn
  ) => ContinuousHistogramBy;

  histogramCategorical: () => CategoricalHistogram;
  histogramCategoricalBy: (by: DataframeColumn) => CategoricalHistogramBy;

  has: (label: LabelType) => boolean;
  ihas: (offset: OffsetType) => boolean;
  indexOf: (value: DataframeValue) => LabelType | undefined;
  iget: (offset: OffsetType) => DataframeValue;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- a legitimate use of any.
export type AnyFunction<T = unknown> = (...args: any[]) => T;

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- a legitimate use of any.
export type HashArgsFunction = (...args: any[]) => string | number;
