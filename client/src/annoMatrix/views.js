/* eslint-disable max-classes-per-file -- Classes are interrelated*/

/*
Views on the annomatrix.  all API here is defined in viewCreators.js and annoMatrix.js.
*/
import clip from "../util/clip";
import AnnoMatrix from "./annoMatrix";
import { _whereCacheCreate } from "./whereCache";
import { _isContinuousType, _getColumnSchema } from "./schema";

class AnnoMatrixView extends AnnoMatrix {
  constructor(viewOf, rowIndex = null) {
    const nObs = rowIndex ? rowIndex.size() : viewOf.nObs;
    super(viewOf.schema, nObs, viewOf.nVar, rowIndex || viewOf.rowIndex);
    this.viewOf = viewOf;
    this.isView = true;
  }

  addObsAnnoCategory(col, category) {
    const newAnnoMatrix = this._clone();
    newAnnoMatrix.viewOf = this.viewOf.addObsAnnoCategory(col, category);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  async removeObsAnnoCategory(col, category, unassignedCategory) {
    const newAnnoMatrix = this._clone();
    newAnnoMatrix.viewOf = await this.viewOf.removeObsAnnoCategory(
      col,
      category,
      unassignedCategory
    );
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  dropObsColumn(col) {
    const newAnnoMatrix = this._clone();
    newAnnoMatrix.viewOf = this.viewOf.dropObsColumn(col);
    newAnnoMatrix._cache.obs = this._cache.obs.dropCol(col);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }
  dropObsmLayout(layout) {
    const newAnnoMatrix = this._clone();    
    try{
      newAnnoMatrix.viewOf = this.viewOf.dropObsmLayout(layout);      
      newAnnoMatrix._cache.emb = this._cache.emb.dropCol(`${layout}_0`);
      newAnnoMatrix._cache.emb = this._cache.emb.dropCol(`${layout}_1`);
    } catch (error){
      console.log(`Ignoring due to: ${error}`)
    }    
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }
  addObsColumn(colSchema, Ctor, value) {
    const newAnnoMatrix = this._clone();
    newAnnoMatrix.viewOf = this.viewOf.addObsColumn(colSchema, Ctor, value);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  renameObsColumn(oldCol, newCol) {
    const newAnnoMatrix = this._clone();
    newAnnoMatrix.viewOf = this.viewOf.renameObsColumn(oldCol, newCol);
    newAnnoMatrix.schema = newAnnoMatrix.viewOf.schema;
    return newAnnoMatrix;
  }

  async setObsColumnValues(col, rowLabels, value) {
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

  async resetObsColumnValues(col, oldValue, newValue) {
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

  addEmbedding(colSchema) {
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
  constructor(viewOf, mapFn) {
    super(viewOf);
    this.mapFn = mapFn;
  }

  async _doLoad(field, query) {
    const df = await this.viewOf._fetch(field, query);
    const dfMapped = df.mapColumns((colData, colIdx) => {
      const colLabel = df.colIndex.getLabel(colIdx);
      const colSchema = _getColumnSchema(this.schema, field, colLabel);
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
  /*
	A view which is a clipped transformation of its parent
	*/
  constructor(viewOf, qmin, qmax) {
    super(viewOf, (field, colLabel, colSchema, colData, df) =>
      _clipAnnoMatrix(field, colLabel, colSchema, colData, df, qmin, qmax)
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
  constructor(viewOf, rowIndex) {
    super(viewOf, rowIndex);
    Object.seal(this);
  }

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
