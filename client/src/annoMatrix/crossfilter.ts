/*
Row crossfilter proxy for an AnnoMatrix.  This wraps Crossfilter,
providing a number of services, and ensuring that the crossfilter and
AnnoMatrix stay in sync:
	- on-demand index creation as data is loaded
	- transparently mapping between queries and crossfilter index names.
  - for mutation of the matrix by user annotations, maintain synchronization
    between Crossfilter and AnnoMatrix.
*/
import Crossfilter from "../util/typedCrossfilter";
import { _getColumnSchema } from "./schema";
import {
  AnnotationColumnSchema,
  Field,
  EmbeddingSchema,
} from "../common/types/schema";
import AnnoMatrix from "./annoMatrix";
import {
  Dataframe,
  DataframeValue,
  DataframeValueArray,
  LabelType,
} from "../util/dataframe";
import { Query } from "./query";
import { AnyArray } from "../common/types/arraytypes";
import { LabelArray } from "../util/dataframe/types";

type ObsDimensionParams =
  | [string, DataframeValueArray, DataframeValueArray]
  | [string, DataframeValueArray]
  | [string, DataframeValueArray, Int32ArrayConstructor]
  | [string, DataframeValueArray, Float32ArrayConstructor];

function _dimensionNameFromDf(field: Field, df: Dataframe): string {
  const colNames = df.colIndex.labels();
  return _dimensionName(field, colNames);
}

function _dimensionName(
  field: Field,
  colNames: LabelType | LabelArray
): string {
  if (!Array.isArray(colNames)) return `${field}/${colNames}`;
  return `${field}/${colNames.join(":")}`;
}

export default class AnnoMatrixObsCrossfilter {
  annoMatrix: AnnoMatrix;

  obsCrossfilter: Crossfilter;

  constructor(
    annoMatrix: AnnoMatrix,
    _obsCrossfilter: Crossfilter | null = null
  ) {
    this.annoMatrix = annoMatrix;
    this.obsCrossfilter =
      _obsCrossfilter || new Crossfilter(annoMatrix._cache.obs);
    this.obsCrossfilter = this.obsCrossfilter.setData(annoMatrix._cache.obs);
  }

  size(): number {
    return this.obsCrossfilter.size();
  }

  /**
  Managing the associated annoMatrix.  These wrappers are necessary to
  make coordinated changes to BOTH the crossfilter and annoMatrix, and
  ensure that all state stays synchronized.

  See API documentation in annoMatrix.js.
  **/
  addObsColumn<T extends DataframeValueArray>(
    colSchema: AnnotationColumnSchema,
    Ctor: new (n: number) => T,
    value: T
  ): AnnoMatrixObsCrossfilter {
    const annoMatrix = this.annoMatrix.addObsColumn(colSchema, Ctor, value);
    const obsCrossfilter = this.obsCrossfilter.setData(annoMatrix._cache.obs);
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  dropObsColumn(col: LabelType): AnnoMatrixObsCrossfilter {
    const annoMatrix = this.annoMatrix.dropObsColumn(col);
    let { obsCrossfilter } = this;
    const dimName = _dimensionName(Field.obs, col);
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  renameObsColumn(
    oldCol: LabelType,
    newCol: LabelType
  ): AnnoMatrixObsCrossfilter {
    const annoMatrix = this.annoMatrix.renameObsColumn(oldCol, newCol);
    const oldDimName = _dimensionName(Field.obs, oldCol);
    const newDimName = _dimensionName(Field.obs, newCol);
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(oldDimName)) {
      obsCrossfilter = obsCrossfilter.renameDimension(oldDimName, newDimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  addObsAnnoCategory(
    col: LabelType,
    category: string
  ): AnnoMatrixObsCrossfilter {
    const annoMatrix = this.annoMatrix.addObsAnnoCategory(col, category);
    const dimName = _dimensionName(Field.obs, col);
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  async removeObsAnnoCategory(
    col: LabelType,
    category: string,
    unassignedCategory: string
  ): Promise<AnnoMatrixObsCrossfilter> {
    const annoMatrix = await this.annoMatrix.removeObsAnnoCategory(
      col,
      category,
      unassignedCategory
    );
    const dimName = _dimensionName(Field.obs, col);
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  async setObsColumnValues(
    col: LabelType,
    rowLabels: Int32Array,
    value: DataframeValue
  ): Promise<AnnoMatrixObsCrossfilter> {
    const annoMatrix = await this.annoMatrix.setObsColumnValues(
      col,
      rowLabels,
      value
    );
    const dimName = _dimensionName(Field.obs, col);
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  async resetObsColumnValues<T extends DataframeValue>(
    col: LabelType,
    oldValue: T,
    newValue: T
  ): Promise<AnnoMatrixObsCrossfilter> {
    const annoMatrix = await this.annoMatrix.resetObsColumnValues(
      col,
      oldValue,
      newValue
    );
    const dimName = _dimensionName(Field.obs, col);
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  addEmbedding(colSchema: EmbeddingSchema): AnnoMatrixObsCrossfilter {
    const annoMatrix = this.annoMatrix.addEmbedding(colSchema);
    return new AnnoMatrixObsCrossfilter(annoMatrix, this.obsCrossfilter);
  }

  /**
   * Drop the crossfilter dimension. Do not change the annoMatrix. Useful when we
   * want to stop tracking the selection state, but aren't sure we want to blow the
   * annomatrix cache.
   */
  dropDimension(field: Field, query: Query): AnnoMatrixObsCrossfilter {
    const { annoMatrix } = this;
    let { obsCrossfilter } = this;
    const keys = annoMatrix
      .getCacheKeys(field, query)
      // @ts-expect-error ts-migrate --- suppressing TS defect (https://github.com/microsoft/TypeScript/issues/44373).
      // Compiler is complaining that expression is not callable on array union types. Remove suppression once fixed.
      .filter((k?: string | number) => k !== undefined);
    const dimName = _dimensionName(field, keys as string[]);
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  /**
  Selection state - API is identical to ImmutableTypedCrossfilter, as these
  are just wrappers to lazy create indices.
  **/

  async select(
    field: Field,
    query: Query,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- TODO revisit: waiting for typings from util/typedCrossfilter
    spec: any
  ): Promise<AnnoMatrixObsCrossfilter> {
    const { annoMatrix } = this;
    let { obsCrossfilter } = this;

    if (!annoMatrix?._cache?.[field]) {
      throw new Error("Unknown field name");
    }
    if (field === "var") {
      throw new Error("unable to obsSelect upon the var dimension");
    }

    // grab the data, so we can grab the index.
    const df = await annoMatrix.fetch(field, query);
    if (!df) {
      throw new Error("Dataframe cannot be `undefined`");
    }
    const dimName = _dimensionNameFromDf(field, df);
    if (!obsCrossfilter.hasDimension(dimName)) {
      // lazy index generation - add dimension when first used
      obsCrossfilter = this._addObsCrossfilterDimension(
        annoMatrix,
        obsCrossfilter,
        field,
        df
      );
    }

    // select
    obsCrossfilter = obsCrossfilter.select(dimName, spec);
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  selectAll(): AnnoMatrixObsCrossfilter {
    /*
		Select all on any dimension in this field.
		*/
    const { annoMatrix } = this;
    const currentDims = this.obsCrossfilter.dimensionNames();
    const obsCrossfilter = currentDims.reduce(
      (xfltr, dim) => xfltr.select(dim, { mode: "all" }),
      this.obsCrossfilter
    );
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  countSelected(): number {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (this.obsCrossfilter.size() === 0) return this.annoMatrix.nObs;
    return this.obsCrossfilter.countSelected();
  }

  allSelectedMask(): Uint8Array {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (
      this.obsCrossfilter.size() === 0 ||
      this.obsCrossfilter.dimensionNames().length === 0
    ) {
      /* fake the mask */
      return new Uint8Array(this.annoMatrix.nObs).fill(1);
    }
    return this.obsCrossfilter.allSelectedMask();
  }

  allSelectedLabels(): LabelArray {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (
      this.obsCrossfilter.size() === 0 ||
      this.obsCrossfilter.dimensionNames().length === 0
    ) {
      return this.annoMatrix.rowIndex.labels();
    }

    const mask = this.obsCrossfilter.allSelectedMask();
    const index = this.annoMatrix.rowIndex.isubsetMask(mask);
    return index.labels();
  }

  fillByIsSelected<T extends AnyArray>(
    array: T,
    selectedValue: T[0],
    deselectedValue: T[0]
  ): T {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (
      this.obsCrossfilter.size() === 0 ||
      this.obsCrossfilter.dimensionNames().length === 0
    ) {
      // @ts-expect-error ts-migrate --- TODO revisit:
      // Type 'Int8Array | Uint8Array | Int16Array | Uint16Array | Int32Array | Uint32Array | Float32Array | Float64Array | unknown[]' is not assignable to type 'T'...
      return array.fill(selectedValue);
    }
    return this.obsCrossfilter.fillByIsSelected(
      array,
      selectedValue,
      deselectedValue
    );
  }

  /**
   ** Private below
   **/

  _addObsCrossfilterDimension(
    annoMatrix: AnnoMatrix,
    obsCrossfilter: Crossfilter,
    field: Field,
    df: Dataframe
  ): Crossfilter {
    if (field === "var") return obsCrossfilter;
    const dimName = _dimensionNameFromDf(field, df);
    const dimParams = this._getObsDimensionParams(field, df);
    obsCrossfilter = obsCrossfilter.setData(annoMatrix._cache.obs);
    // @ts-expect-error ts-migrate --- TODO revisit:
    // `...dimParams`: A spread argument must either have a tuple type or be passed to a rest parameter.
    obsCrossfilter = obsCrossfilter.addDimension(dimName, ...dimParams);
    return obsCrossfilter;
  }

  _getColumnBaseType(field: Field, col: LabelType): string {
    /* Look up the primitive type for this field/col */
    const colSchema = _getColumnSchema(this.annoMatrix.schema, field, col);
    return colSchema.type;
  }

  _getObsDimensionParams(
    field: Field,
    df: Dataframe
  ): ObsDimensionParams | undefined {
    /* return the crossfilter dimensiontype type and params for this field/dataframe */

    if (field === Field.emb) {
      /* assumed to be 2D */
      return ["spatial", df.icol(0).asArray(), df.icol(1).asArray()];
    }

    /* assumed to be 1D */
    const col = df.icol(0);
    const colName = df.colIndex.getLabel(0);
    // @ts-expect-error --- TODO revisit:
    // `colName` Argument of type 'LabelType | undefined' is not assignable to parameter of type 'LabelType'. Type 'undefined' is not assignable to type 'LabelType'.
    const type = this._getColumnBaseType(field, colName);
    if (type === "string" || type === "categorical" || type === "boolean") {
      return ["enum", col.asArray()];
    }
    if (type === "int32") {
      return ["scalar", col.asArray(), Int32Array];
    }
    if (type === "float32") {
      return ["scalar", col.asArray(), Float32Array];
    }
    // Currently not supporting boolean and categorical types.
    console.error(
      `Warning - unknown metadata schema (${type}) for field ${field} ${colName}.`
    );
    // skip it - we don't know what to do with this type

    return undefined;
  }
}
