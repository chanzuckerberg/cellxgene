/* eslint-disable max-classes-per-file -- Classes are interrelated*/

/*
Views on the annomatrix.  all API here is defined in viewCreators.js and annoMatrix.js.
*/
import clip from "../util/clip";
import AnnoMatrix from "./annoMatrix";
import { _whereCacheCreate } from "./whereCache";
import { _isContinuousType, _getColumnSchema } from "./schema";

class AnnoMatrixView extends AnnoMatrix {
  _cache: any;

  constructor(viewOf: any, rowIndex = null) {
    // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
    const nObs = rowIndex ? rowIndex.size() : viewOf.nObs;
    super(viewOf.schema, nObs, viewOf.nVar, rowIndex || viewOf.rowIndex);
    this.viewOf = viewOf;
    this.isView = true;
  }

  addObsAnnoCategory(col: any, category: any) {
    const newAnnoMatrix = this._clone();
    newAnnoMatrix.viewOf = this.viewOf.addObsAnnoCategory(col, category);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  async removeObsAnnoCategory(
    col: any,
    category: any,
    unassignedCategory: any
  ) {
    const newAnnoMatrix = this._clone();
    newAnnoMatrix.viewOf = await this.viewOf.removeObsAnnoCategory(
      col,
      category,
      unassignedCategory
    );
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  dropObsColumn(col: any) {
    const newAnnoMatrix = this._clone();
    newAnnoMatrix.viewOf = this.viewOf.dropObsColumn(col);
    newAnnoMatrix._cache.obs = this._cache.obs.dropCol(col);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  addObsColumn(colSchema: any, Ctor: any, value: any) {
    const newAnnoMatrix = this._clone();
    newAnnoMatrix.viewOf = this.viewOf.addObsColumn(colSchema, Ctor, value);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  renameObsColumn(oldCol: any, newCol: any) {
    const newAnnoMatrix = this._clone();
    newAnnoMatrix.viewOf = this.viewOf.renameObsColumn(oldCol, newCol);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  async setObsColumnValues(col: any, rowLabels: any, value: any) {
    const newAnnoMatrix = this._clone();
    newAnnoMatrix.viewOf = await this.viewOf.setObsColumnValues(
      col,
      rowLabels,
      value
    );
    newAnnoMatrix._cache.obs = this._cache.obs.dropCol(col);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  async resetObsColumnValues(col: any, oldValue: any, newValue: any) {
    const newAnnoMatrix = this._clone();
    newAnnoMatrix.viewOf = await this.viewOf.resetObsColumnValues(
      col,
      oldValue,
      newValue
    );
    newAnnoMatrix._cache.obs = this._cache.obs.dropCol(col);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  addEmbedding(colSchema: any) {
    const newAnnoMatrix = this._clone();
    newAnnoMatrix.viewOf = this.viewOf.addEmbedding(colSchema);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }
}

class AnnoMatrixMapView extends AnnoMatrixView {
  /*
	A view which knows how to transform its data.
	*/
  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'viewOf' implicitly has an 'any' type.
  constructor(viewOf, mapFn) {
    super(viewOf);
    (this as any).mapFn = mapFn;
  }

  // @ts-expect-error ts-migrate(2416) FIXME: Property '_doLoad' in type 'AnnoMatrixMapView' is ... Remove this comment to see the full error message
  async _doLoad(field, query) {
    const df = await this.viewOf._fetch(field, query);
    // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'colData' implicitly has an 'any' type.
    const dfMapped = df.mapColumns((colData, colIdx) => {
      const colLabel = df.colIndex.getLabel(colIdx);
      const colSchema = _getColumnSchema(this.schema, field, colLabel);
      return (this as any).mapFn(field, colLabel, colSchema, colData, df);
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
  /*
	A view which is a clipped transformation of its parent
	*/
  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'viewOf' implicitly has an 'any' type.
  constructor(viewOf, qmin, qmax) {
    // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'field' implicitly has an 'any' type.
    super(viewOf, (field, colLabel, colSchema, colData, df) =>
      _clipAnnoMatrix(field, colLabel, colSchema, colData, df, qmin, qmax)
    );
    (this as any).isClipped = true;
    (this as any).clipRange = [qmin, qmax];
    Object.seal(this);
  }
}

export class AnnoMatrixRowSubsetView extends AnnoMatrixView {
  /*
	A view which is a subset of total rows.
	*/
  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'viewOf' implicitly has an 'any' type.
  constructor(viewOf, rowIndex) {
    super(viewOf, rowIndex);
    Object.seal(this);
  }

  // @ts-expect-error ts-migrate(2416) FIXME: Property '_doLoad' in type 'AnnoMatrixRowSubsetVie... Remove this comment to see the full error message
  async _doLoad(field, query) {
    const df = await this.viewOf._fetch(field, query);

    // don't try to row-subset the var dimension.
    if (field === "var") {
      return [null, df];
    }

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

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'field' implicitly has an 'any' type.
function _clipAnnoMatrix(field, colLabel, colSchema, colData, df, qmin, qmax) {
  /* only clip obs and var scalar columns */
  if (field !== "obs" && field !== "X") return colData;
  if (!_isContinuousType(colSchema)) return colData;
  if (qmin < 0) qmin = 0;
  if (qmax > 1) qmax = 1;
  if (qmin === 0 && qmax === 1) return colData;

  const quantiles = df.col(colLabel).summarize().percentiles;
  const lower = quantiles[100 * qmin];
  const upper = quantiles[100 * qmax];
  const clippedData = clip(colData.slice(), lower, upper, Number.NaN);
  return clippedData;
}

/* eslint-enable max-classes-per-file -- enable*/
