/* eslint-disable max-classes-per-file -- Classes are interrelated*/

import { doBinaryRequest, fetchResult, dubEncURIComp } from "./fetchHelpers";
import { matrixFBSToDataframe } from "../stateManager/matrix";
import {
  addObsAnnoColumn,
  removeObsAnnoColumn,
  addObsAnnoCategory,
  removeObsAnnoCategory,
} from "../stateManager/schemaHelpers";
import { Dataframe, IdentityInt32Index } from "../dataframe";
import {
  getColumnDimensionNames,
  getColumnSchema,
  schemaColumns,
  indexEntireSchema,
  isContinuousType,
} from "./schema";
import { whereCacheCreate, whereCacheGet, whereCacheMerge } from "./whereCache";
import clip from "../clip";
import { isArrayOrTypedArray } from "../typeHelpers";

class AnnoMatrix {
  static fields() {
    return ["obs", "var", "emb", "X"];
  }

  constructor(schema, nObs, nVar, rowIndex = null) {
    this.schema = indexEntireSchema(schema);
    this.nObs = nObs;
    this.nVar = nVar;
    this.rowIndex = rowIndex || new IdentityInt32Index(nObs);
    this.isView = false;
    this.viewOf = undefined;

    /*
		Private:

		These are caches - lazily loaded. The only guarantee is that if they
		are loaded, they will conform to the schema & dimensionality constraints.

		Do NOT use directly - instead, use the fetch() and preload() API.
		*/
    this.X = Dataframe.empty(this.rowIndex);
    this.obs = Dataframe.empty(this.rowIndex);
    this.var = Dataframe.empty(this.rowIndex);
    this.emb = Dataframe.empty(this.rowIndex);

    /* private */
    this._pendingLoad = {
      obs: {},
      var: {},
      emb: {},
      X: {},
    };
    this._whereCache = {};
  }

  /**
   ** Schema helper/accessors
   **/
  getMatrixColumns(field) {
    return schemaColumns(this.schema, field);
  }

  static getMatrixFields() {
    return AnnoMatrix.fields();
  }

  getColumnSchema(field, col) {
    return getColumnSchema(this.schema, field, col);
  }

  getColumnDimensions(field, col) {
    return getColumnDimensionNames(this.schema, field, col);
  }

  /**
   ** Load / read interfaces
   **/
  fetch(field, q) {
    /*
		Return the given query on a single matrix field as a single dataframe.
		Currently supports ONLY full column query.

		Returns a Promise for the dataframe.

		Field must be one of the matrix fields: 'obs', 'var', 'X', 'emb'.  Value
		represents the underlying object upon which the query is occuring.

		Query is one of:
			* a string, representing a single column name from the field, eg,
				"n_genes"
			* an object, containing an advanced query
			* an array, containing one or more of the above.

		Columns may have more than one dimension, and all will be fetched
		and returned together.  This is most common in an embedding.

		A value query allows for fetching based upon the value in another 
		field/column, _on the same dimension_ (currently only on the var 
		dimension, as only full column fetches are supported). The query
		is a single value filter:
			{ "field name": [
				{name: "column name", values: [ list of values ]}
			]}
		One and only one value filter is allowed.

		Examples:
			fetch("obs", "n_genes") // get the n_genes column
			fetch("var", this.schema.annotations.var.index)	// get the var index
			fetch("obs", ["n_genes", "louvain"]) // multiple columns
			fetch("X", { 
				where: {field: "var", column: this.schema.annotations.var.index, value: "TYMP"}
			})

		The value query is a recodification and subset of the server REST API 
		value filter JSON.  Range queries and multiple filters are not currently
		supported.

		Returns a Promise for a dataframe containing the requested columns.
		*/
    return fetchResult(this._fetch(field, q));
  }

  prefetch(field, q) {
    /*
		Start a data fetch & cache fill.  Identical to fetch() except it does
		not return a value.
		*/
    this._fetch(field, q);
    return undefined;
  }

  /**
   ** View creation interfaces:
   **/
  isubsetMask(rowMask) {
    /*
		Subset on row based upon mask
		*/
    return this.isubset(maskToList(rowMask));
  }

  isubset(rowOffsets) {
    /*
		subset based on offset
		*/
    const rowIndex = this.rowIndex.isubset(rowOffsets);
    return new AnnoMatrixRowSubsetView(this, rowIndex);
  }

  subset(rowLabels) {
    /*
		subset based on labels
		*/
    const rowIndex = this.rowIndex.subset(rowLabels);
    return new AnnoMatrixRowSubsetView(this, rowIndex);
  }

  clip(qmin, qmax) {
    /*
		Create a view that clips all continuous data to the [min, max] range.
		The matrix shape does not change, but the continuous values outside the
		specified range will become a NaN.
		*/
    return new AnnoMatrixClipView(this, qmin, qmax);
  }

  /**
   ** private below
   **/
  async _fetch(field, q) {
    if (!AnnoMatrix.fields().includes(field)) return undefined;
    const queries = Array.isArray(q) ? q : [q];

    /* find any query not already cached */
    const uncachedQueries = queries.filter((query) =>
      whereCacheGet(this._whereCache, this.schema, field, query).some(
        (cacheKey) => cacheKey === undefined || !this[field].hasCol(cacheKey)
      )
    );

    /* load uncached queries */
    if (uncachedQueries.length > 0) {
      await Promise.all(
        uncachedQueries.map((query) =>
          this._getPendingLoad(field, query, async (_field, _query) => {
            /* fetch, then index */
            const [whereCacheUpdate, df] = await this._doLoad(_field, _query);
            this[_field] = this[_field].withColsFrom(df);
            this._whereCache = whereCacheMerge(
              this._whereCache,
              whereCacheUpdate
            );
          })
        )
      );
    }

    /* everything we need is in the cache, so just cherry-pick requested columns */
    const requestedCacheKeys = queries
      .map((query) =>
        whereCacheGet(this._whereCache, this.schema, field, query)
      )
      .flat();
    return this[field].subset(null, requestedCacheKeys);
  }

  async _getPendingLoad(field, query, fetchFn) {
    const key = queryCacheKey(field, query);
    if (!this._pendingLoad[field][key]) {
      this._pendingLoad[field][key] = fetchFn(field, query);
      await this._pendingLoad[field][key];
      delete this._pendingLoad[field][key];
    }
    return this._pendingLoad[field][key];
  }
}

/*
Base matrix, which proxies for server.
*/
export class AnnoMatrixLoader extends AnnoMatrix {
  constructor(baseURL, schema) {
    const { nObs, nVar } = schema.dataframe;
    super(schema, nObs, nVar);
    this.baseURL = baseURL;
  }

  async _doLoad(field, query) {
    /*
		_doLoad - evaluates the query against the field. Returns:
			* whereCache update: column query map mapping the query to the column labels
			* Dataframe contianing the new colums (one per dimension)
		*/
    let urlQuery;
    let urlBase;

    switch (field) {
      case "obs":
      case "var": {
        urlBase = `${this.baseURL}/annotations/${field}`;
        urlQuery = encodeQuery("annotation-name", query);
        break;
      }
      case "X": {
        urlBase = `${this.baseURL}/data/var`;
        urlQuery = encodeQuery(undefined, query);
        break;
      }
      case "emb": {
        urlBase = `${this.baseURL}/layout/obs`;
        urlQuery = encodeQuery("layout-name", query);
        break;
      }
      default:
        throw new Error("Unknown field name");
    }

    const url = `${urlBase}?${urlQuery}`;
    const buffer = await doBinaryRequest(url);
    const result = matrixFBSToDataframe(buffer);
    if (!result || result.isEmpty()) throw Error("Unknown field/col");

    const whereCacheUpdate = whereCacheCreate(
      field,
      query,
      result.colIndex.labels()
    );
    return [whereCacheUpdate, result];
  }

  addObsAnnoCategory(col, category) {
    const o = clone(this);
    o.schema = addObsAnnoCategory(this.schema, col, category);
    return o;
  }

  removeObsAnnoCategory(col, category) {
    const o = clone(this);
    o.schema = removeObsAnnoCategory(this.schema, col, category);
    return o;
  }

  dropObsColumn(col) {
    /*
		drop column from field
		*/
    const field = "obs";
    const colSchema = getColumnSchema(this.schema, field, col);
    if (!colSchema.writable) {
      throw new Error("Not a writable obs column");
    }

    const o = clone(this);
    o.obs = this.obs.dropCol(col);
    o.schema = removeObsAnnoColumn(this.schema, col);
    return o;
  }

  addObsColumn(colSchema, Ctor, value) {
    /*
		add a column to field, initializing with value.  Value may 
    be one of:
      * an array of values
      * a primitive type, including null or undefined.
    If an array, it must be of same size as nObs and same type as Ctor
		*/
    const field = "obs";
    const col = colSchema.name;
    if (getColumnSchema(this.schema, field, col) || this[field].hasCol(col)) {
      throw new Error("column already exists");
    }

    const o = clone(this);
    let data;
    if (isArrayOrTypedArray(value)) {
      if (value.constructor !== Ctor)
        throw new Error("Mismatched value array type");
      if (value.length !== this.nObs)
        throw new Error("Value array has incorrect length");
      data = value.slice();
    } else {
      data = new Ctor(this.nObs).fill(value);
    }
    o.obs = this.obs.withCol(col, data);
    o.schema = addObsAnnoColumn(this.schema, col, {
      ...colSchema,
      writable: true,
    });
    return o;
  }

  renameObsColumn(oldCol, newCol) {
    /*
    Rename the obs oldColName to newColName.  oldCol must be writable.
    */
    const field = "obs";
    const oldColSchema = getColumnSchema(this.schema, field, oldCol);
    if (!oldColSchema.writable) {
      throw new Error("Not a writable obs column");
    }

    const value = this[field].hasCol(oldCol)
      ? this[field].col(oldCol).asArray()
      : undefined;
    return this.dropObsColumn(oldCol).addObsColumn(
      {
        ...oldColSchema,
        name: newCol,
      },
      value.constructor,
      value
    );
  }

  setObsColumnValues(col, rowLabels, value) {
    /*
		Set all rows in rowLabels to value.
		*/
    const field = "obs";
    const colSchema = getColumnSchema(this.schema, field, col);
    if (!colSchema.writable) {
      throw new Error("Not a writable obs column");
    }

    const rowIndices = this.rowIndex.getOffsets(rowLabels);
    const data = this.obs.col(col).asArray().slice();
    for (let i = 0, len = rowIndices.length; i < len; i += 1) {
      data[rowIndices[i]] = value;
    }

    const o = clone(this);
    o.obs = this.obs.replaceColData(col, data);
    if (colSchema.type === "categorical") {
      o.schema = addObsAnnoCategory(this.schema, col, value);
    }
    return o;
  }
}

class AnnoMatrixView extends AnnoMatrix {
  constructor(viewOf, rowIndex = null) {
    const nObs = rowIndex ? rowIndex.size() : viewOf.nObs;
    super(viewOf.schema, nObs, viewOf.nVar, rowIndex || viewOf.rowIndex);
    this.viewOf = viewOf;
    this.isView = true;
  }

  addObsAnnoCategory(col, category) {
    const o = clone(this);
    o.viewOf = this.viewOf.addObsAnnoCategory(col, category);
    o.schema = o.viewOf.schema;
    return o;
  }

  removeObsAnnoCategory(col, category) {
    const o = clone(this);
    o.viewOf = this.viewOf.removeObsAnnoCategory(col, category);
    o.schema = o.viewOf.schema;
    return o;
  }

  dropObsColumn(col) {
    /*
		drop column from field
		*/
    // const field = "obs";
    // const colSchema = getColumnSchema(this.schema, field, col);
    // if (!colSchema.writable) {
    //   throw new Error("Not a writable obs column");
    // }

    const o = clone(this);
    o.viewOf = this.viewOf.dropObsColumn(col);
    o.obs = this.obs.dropCol(col);
    o.schema = o.viewOf.schema;
    return o;
  }

  addObsColumn(colSchema, ctor, value) {
    /*
		add a column to field, initializing with value
		*/
    // const field = "obs";
    // const col = colSchema.name;
    // if (getColumnSchema(this.schema, field, col) || this[field].hasCol(col)) {
    //   throw new Error("column already exists");
    // }

    const o = clone(this);
    o.viewOf = this.viewOf.addObsColumn(colSchema, ctor, value);
    o.schema = o.viewOf.schema;
    return o;
  }

  setObsColumnValues(col, rowLabels, value) {
    /*
		set values
		*/
    // const field = "obs";
    // const colSchema = getColumnSchema(this.schema, field, col);
    // if (!colSchema.writable) {
    //   throw new Error("Not a writable obs column");
    // }

    const o = clone(this);
    o.viewOf = this.viewOf.setObsColumnValues(col, rowLabels, value);
    o.obs = this.obs.dropCol(col);
    o.schema = o.viewOf.schema;
    return o;
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
      const colSchema = getColumnSchema(this.schema, field, colLabel);
      return this.mapFn(field, colLabel, colSchema, colData, df);
    });
    const whereCacheUpdate = whereCacheCreate(
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
      clipAnnoMatrix(field, colLabel, colSchema, colData, df, qmin, qmax)
    );
    this.isClipped = true;
    this.clipRange = [qmin, qmax];
  }
}

export class AnnoMatrixRowSubsetView extends AnnoMatrixView {
  /*
	A view which is a subset of total rows.
	*/
  async _doLoad(field, query) {
    const df = await this.viewOf._fetch(field, query);

    // don't try to row-subset the var dimension.
    if (field === "var") {
      return [null, df];
    }

    const dfSubset = df.subset(null, null, this.rowIndex);
    const whereCacheUpdate = whereCacheCreate(
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

function clipAnnoMatrix(field, colLabel, colSchema, colData, df, qmin, qmax) {
  /* only clip obs and var scalar columns */
  if (field !== "obs" && field !== "X") return colData;
  if (!isContinuousType(colSchema)) return colData;
  if (qmin < 0) qmin = 0;
  if (qmax > 1) qmax = 1;
  if (qmin === 0 && qmax === 1) return colData;

  const quantiles = df.col(colLabel).summarize().percentiles;
  const lower = quantiles[100 * qmin];
  const upper = quantiles[100 * qmax];
  const clippedData = clip(colData.slice(), lower, upper, Number.NaN);
  return clippedData;
}

/* convert masks to lists - method wastes space, but is fast */
function maskToList(mask) {
  if (!mask) {
    return null;
  }
  const list = new Int32Array(mask.length);
  let elems = 0;
  for (let i = 0, l = mask.length; i < l; i += 1) {
    if (mask[i]) {
      list[elems] = i;
      elems += 1;
    }
  }
  return new Int32Array(list.buffer, 0, elems);
}

function clone(orig) {
  return Object.assign(Object.create(Object.getPrototypeOf(orig)), orig);
}

function queryCacheKey(field, query) {
  if (typeof query === "object") {
    const { field: queryField, column: queryColumn, value: queryValue } = query;
    return `${field}/${queryField}/${queryColumn}/${queryValue}`;
  }
  return `${field}/${query}`;
}

function encodeQuery(colKey, q) {
  if (typeof q === "object") {
    const { field: queryField, column: queryColumn, value: queryValue } = q;
    return `${dubEncURIComp(queryField)}:${dubEncURIComp(
      queryColumn
    )}=${dubEncURIComp(queryValue)}`;
  }
  if (!colKey) throw new Error("Unsupported query by name");
  return `${colKey}=${encodeURIComponent(q)}`;
}
/* eslint-enable max-classes-per-file -- enable*/
