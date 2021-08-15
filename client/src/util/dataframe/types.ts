import { TypedArray } from "../../common/types/arraytypes";

export type LabelType = number | string;

type CommonProps<A, B> = {
  [K in keyof A & keyof B]: A[K] | B[K];
};
export type GenericLabelArray<T> = CommonProps<Array<T>, Int32Array>;
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

export type DataframeValueArray = DataframeValue[] | TypedArray;

export type DataframeColumnGetter = (
  label: LabelType
) => DataframeValue | undefined;

/**
 * Interface representing a Dataframe column.  Eg, returned by
 * Dataframe.col().
 */
export interface DataframeColumn extends DataframeColumnGetter {
  /**
   * __id is unique per Dataframe and DataframeColumn, and is used as a memoization key.
   */
  readonly __id: string;

  /**
   * Boolean indicating if the underlying data supports continuous operations, eg,
   * summarizeContinuous.
   */
  isContinuous: boolean;

  /**
   * Return underlying column data as an array-like object.
   */
  asArray: () => DataframeValueArray;

  /**
   * Continuous data summary.  Will throw if !isContinuous.
   */
  summarizeContinuous: () => ContinuousColumnSummary;

  /**
   * Categorical data summary.
   */
  summarizeCategorical: () => CategoricalColumnSummary;

  /**
   * Continuous bin/histogram.  Will throw if !isContinuous.
   * @param bins - array of bin boundary fractions, in range [0., 1.]
   * @param domain - data domain [min, max]
   */
  histogramContinuous: (
    bins: number,
    domain: [number, number]
  ) => ContinuousHistogram;

  /**
   * Continuous bin/histogram, grouped by another categorical column.  Will throw if !isContinuous.
   * @param bins - array of bin boundary fractions, in range [0., 1.]
   * @param domain - data domain [min, max]
   * @param by - group by categorical column
   */
  histogramContinuousBy: (
    bins: number,
    domain: [number, number],
    by: DataframeColumn
  ) => ContinuousHistogramBy;

  /**
   * Categorical bin/histogram.
   */
  histogramCategorical: () => CategoricalHistogram;

  /**
   * Categorical bin/histogram grouped by another column.
   * @param by - group by categorical column
   */
  histogramCategoricalBy: (by: DataframeColumn) => CategoricalHistogramBy;

  /**
   * Return true if the column contains the row label.
   */
  has: (rlabel: LabelType) => boolean;

  /**
   * Return true if the column contains the row offset. Identical to
   * (offset >= 0 && offset < dataframe.length)
   */
  ihas: (offset: OffsetType) => boolean;

  /**
   * Return index of the value, as a _label_. Returns undefined if
   * not present.  *NOTE*: unlike Array.indexOf, does not return an
   * offset.
   */
  indexOf: (value: DataframeValue) => LabelType | undefined;

  /**
   * Return the value at the given offset, or undefined if not present.
   */
  iget: (offset: OffsetType) => DataframeValue | undefined;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- a legitimate use of any.
export type AnyFunction<T = unknown> = (...args: any[]) => T;

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- a legitimate use of any.
export type HashArgsFunction = (...args: any[]) => string | number;
