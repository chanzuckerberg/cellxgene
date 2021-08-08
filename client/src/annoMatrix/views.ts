/* eslint-disable max-classes-per-file -- Classes are interrelated*/

/*
Views on the annomatrix.  all API here is defined in viewCreators.js and annoMatrix.js.
*/
import clip from "../util/clip";
import AnnoMatrix, { ObsColumnValue } from "./annoMatrix";
import { _whereCacheCreate, WhereCache } from "./whereCache";
import { _isContinuousType, _getColumnSchema } from "./schema";
import {
  Dataframe,
  DenseInt32Index,
  IdentityInt32Index,
  KeyIndex,
} from "../util/dataframe";
import { Query } from "./query";
import {
  AnnotationColumnSchema,
  ArraySchema,
  Field,
  EmbeddingSchema,
} from "../common/types/schema";
import { AnyArray } from "../common/types/arraytypes";

type MapFn = (
  field: Field,
  colLabel: string,
  colSchema: ArraySchema,
  colData: string[],
  df: Dataframe
) => string[];

abstract class AnnoMatrixView extends AnnoMatrix {
  constructor(
    viewOf: AnnoMatrix,
    rowIndex: DenseInt32Index | IdentityInt32Index | KeyIndex | null = null
  ) {
    const nObs = rowIndex ? rowIndex.size() : viewOf.nObs;
    super(viewOf.schema, nObs, viewOf.nVar, rowIndex || viewOf.rowIndex);
    this.viewOf = viewOf;
    this.isView = true;
  }

  addObsAnnoCategory(col: string, category: string): AnnoMatrix {
    const newAnnoMatrix = this._clone();
    // @ts-expect-error ts-migrate --- TODO revisit:
    // `this.viewOf`: Object is possibly 'undefined'.
    newAnnoMatrix.viewOf = this.viewOf.addObsAnnoCategory(col, category);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  async removeObsAnnoCategory(
    col: string,
    category: string,
    unassignedCategory: string
  ): Promise<AnnoMatrix> {
    const newAnnoMatrix = this._clone();
    // @ts-expect-error ts-migrate --- TODO revisit: see comments on L46
    newAnnoMatrix.viewOf = await this.viewOf.removeObsAnnoCategory(
      col,
      category,
      unassignedCategory
    );
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  dropObsColumn(col: string): AnnoMatrix {
    const newAnnoMatrix = this._clone();
    // @ts-expect-error ts-migrate --- TODO revisit: see comments on L46
    newAnnoMatrix.viewOf = this.viewOf.dropObsColumn(col);
    newAnnoMatrix._cache.obs = this._cache.obs.dropCol(col);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  addObsColumn<T extends ObsColumnValue>(
    colSchema: AnnotationColumnSchema,
    Ctor: new (n: number) => T[],
    value: T | T[]
  ): AnnoMatrix {
    const newAnnoMatrix = this._clone();
    // @ts-expect-error ts-migrate --- TODO revisit: see comments on L46
    newAnnoMatrix.viewOf = this.viewOf.addObsColumn(colSchema, Ctor, value);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  renameObsColumn(oldCol: string, newCol: string): AnnoMatrix {
    const newAnnoMatrix = this._clone();
    // @ts-expect-error ts-migrate --- TODO revisit: see comments on L46
    newAnnoMatrix.viewOf = this.viewOf.renameObsColumn(oldCol, newCol);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  async setObsColumnValues(
    col: string,
    rowLabels: Int32Array,
    value: ObsColumnValue
  ): Promise<AnnoMatrix> {
    const newAnnoMatrix = this._clone();
    // @ts-expect-error ts-migrate --- TODO revisit: see comments on L46
    newAnnoMatrix.viewOf = await this.viewOf.setObsColumnValues(
      col,
      rowLabels,
      value
    );
    newAnnoMatrix._cache.obs = this._cache.obs.dropCol(col);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  async resetObsColumnValues<T extends ObsColumnValue>(
    col: string,
    oldValue: T,
    newValue: T
  ): Promise<AnnoMatrix> {
    const newAnnoMatrix = this._clone();
    // @ts-expect-error ts-migrate --- TODO revisit: see comments on L46
    newAnnoMatrix.viewOf = await this.viewOf.resetObsColumnValues(
      col,
      oldValue,
      newValue
    );
    newAnnoMatrix._cache.obs = this._cache.obs.dropCol(col);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  addEmbedding(colSchema: EmbeddingSchema): AnnoMatrix {
    const newAnnoMatrix = this._clone();
    // @ts-expect-error ts-migrate --- TODO revisit: see comments on L46
    newAnnoMatrix.viewOf = this.viewOf.addEmbedding(colSchema);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }
}

class AnnoMatrixMapView extends AnnoMatrixView {
  mapFn: MapFn;

  /*
    A view which knows how to transform its data.
    */
  constructor(viewOf: AnnoMatrix, mapFn: MapFn) {
    super(viewOf);
    this.mapFn = mapFn;
  }

  async _doLoad(
    field: Field,
    query: Query
  ): Promise<[WhereCache | null, Dataframe]> {
    // @ts-expect-error ts-migrate --- TODO revisit:
    // `this.viewOf`: Object is possibly 'undefined'.
    const df = await this.viewOf._fetch(field, query);
    // @ts-expect-error ts-migrate --- TODO revisit:
    // `df`: Object is possibly 'undefined'.
    const dfMapped = df.mapColumns((colData: string[], colIdx: number) => {
      // @ts-expect-error ts-migrate --- TODO revisit:
      // `df`: Object is possibly 'undefined'.
      const colLabel = df.colIndex.getLabel(colIdx);
      const colSchema = _getColumnSchema(this.schema, field, colLabel);
      // @ts-expect-error ts-migrate --- TODO revisit:
      // `df`: Object is possibly 'undefined'.
      return this.mapFn(field, colLabel, colSchema, colData, df);
    });
    const whereCacheUpdate = _whereCacheCreate(
      field,
      query,
      dfMapped.colIndex.labels()
    );
    return [whereCacheUpdate, dfMapped];
  }
}

export class AnnoMatrixClipView extends AnnoMatrixMapView {
  clipRange: [number, number];

  isClipped: boolean;

  /*
    A view which is a clipped transformation of its parent
    */
  constructor(viewOf: AnnoMatrix, qmin: number, qmax: number) {
    super(
      viewOf,
      (
        field: Field,
        colLabel: string,
        colSchema: ArraySchema,
        colData: string[],
        df: Dataframe
      ) => _clipAnnoMatrix(field, colLabel, colSchema, colData, df, qmin, qmax)
    );
    this.isClipped = true;
    this.clipRange = [qmin, qmax];
    Object.seal(this);
  }
}

export class AnnoMatrixRowSubsetView extends AnnoMatrixView {
  /*
    A view which is a subset of total rows.
    */
  constructor(
    viewOf: AnnoMatrix,
    rowIndex: DenseInt32Index | IdentityInt32Index | KeyIndex
  ) {
    super(viewOf, rowIndex);
    Object.seal(this);
  }

  async _doLoad(
    field: Field,
    query: Query
  ): Promise<[WhereCache | null, Dataframe]> {
    // @ts-expect-error ts-migrate --- TODO revisit:
    // `this.viewOf`: Object is possibly 'undefined'.
    const df = await this.viewOf._fetch(field, query);

    // don't try to row-subset the var dimension.
    if (field === Field.var) {
      // @ts-expect-error ts-migrate --- TODO revisit:
      // `df`:  Type 'Dataframe | undefined' is not assignable to type 'Dataframe'. Type 'undefined' is not assignable to type 'Dataframe'.
      return [null, df];
    }

    // @ts-expect-error ts-migrate --- TODO revisit:
    // `df`: Object is possibly 'undefined'.
    const dfSubset = df.subset(null, null, this.rowIndex);
    const whereCacheUpdate = _whereCacheCreate(
      field,
      query,
      dfSubset.colIndex.labels()
    );
    return [whereCacheUpdate, dfSubset];
  }
}

/*
Utility functions below
*/

function _clipAnnoMatrix<T extends AnyArray>(
  field: Field,
  colLabel: string,
  colSchema: ArraySchema,
  colData: T,
  df: Dataframe,
  qmin: number,
  qmax: number
): T {
  /* only clip obs and var scalar columns */
  if (field !== Field.obs && field !== Field.X) return colData;
  if (!_isContinuousType(colSchema)) return colData;
  if (qmin < 0) qmin = 0;
  if (qmax > 1) qmax = 1;
  if (qmin === 0 && qmax === 1) return colData;

  const quantiles = df.col(colLabel).summarizeContinuous().percentiles;
  const lower = quantiles[100 * qmin];
  const upper = quantiles[100 * qmax];
  const clippedData = clip(colData.slice(), lower, upper, Number.NaN);
  return clippedData;
}

/* eslint-enable max-classes-per-file -- enable*/
