import { Dataframe, IdentityInt32Index } from "../util/dataframe";
import {
  _getColumnDimensionNames,
  _getColumnSchema,
  _schemaColumns,
  _getWritableColumns,
} from "./schema";
import { _fetchResult } from "./fetchHelpers";
import { indexEntireSchema } from "../util/stateManager/schemaHelpers";
import { _whereCacheGet, _whereCacheMerge } from "./whereCache";

export default class AnnoMatrix {
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
    this._gcInfo = {
      obs: new Map(),
      var: new Map(),
      emb: new Map(),
      X: new Map(),
    };
  }

  /**
   ** Schema helper/accessors
   **/
  getMatrixColumns(field) {
    return _schemaColumns(this.schema, field);
  }

  // eslint-disable-next-line class-methods-use-this -- need to be able to call this on instances
  getMatrixFields() {
    return AnnoMatrix.fields();
  }

  getColumnSchema(field, col) {
    return _getColumnSchema(this.schema, field, col);
  }

  getColumnDimensions(field, col) {
    return _getColumnDimensionNames(this.schema, field, col);
  }

  /**
   ** General utility methods
   **/
  base() {
    /*
    return the base of view
    */
    let annoMatrix = this;
    while (annoMatrix.isView) annoMatrix = annoMatrix.viewOf;
    return annoMatrix;
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
    return _fetchResult(this._fetch(field, q));
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
   ** private below
   **/
  _resolveCachedQueries(field, queries) {
    return queries
      .map((query) =>
        _whereCacheGet(this._whereCache, this.schema, field, query).filter(
          (cacheKey) => cacheKey !== undefined && this[field].hasCol(cacheKey)
        )
      )
      .flat();
  }

  async _fetch(field, q) {
    if (!AnnoMatrix.fields().includes(field)) return undefined;
    const queries = Array.isArray(q) ? q : [q];

    /* find cached columns we need, and GC the rest */
    const cachedColumns = this._resolveCachedQueries(field, queries);
    this._gcFetchCleanup(field, cachedColumns);

    /* find any query not already cached */
    const uncachedQueries = queries.filter((query) =>
      _whereCacheGet(this._whereCache, this.schema, field, query).some(
        (cacheKey) => cacheKey === undefined || !this[field].hasCol(cacheKey)
      )
    );

    /* load uncached queries */
    if (uncachedQueries.length > 0) {
      await Promise.all(
        uncachedQueries.map((query) =>
          this._getPendingLoad(field, query, async (_field, _query) => {
            /* fetch, then index.  _doLoad is subclass interface */
            const [whereCacheUpdate, df] = await this._doLoad(_field, _query);
            this[_field] = this[_field].withColsFrom(df);
            this._whereCache = _whereCacheMerge(
              this._whereCache,
              whereCacheUpdate
            );
          })
        )
      );
    }

    /* everything we need is in the cache, so just cherry-pick requested columns */
    const requestedCacheKeys = this._resolveCachedQueries(field, queries);
    const response = this[field].subset(null, requestedCacheKeys);
    this._gcUpdateStats(field, response);
    return response;
  }

  async _getPendingLoad(field, query, fetchFn) {
    /*
    Given a query on a field, ensure that we only have a single outstanding
    fetch at any given time.  If multiple requests occur while a fetch is
    outstanding, just wait for the original.

    This is implemented by returning a promise that will await the singular
    fetch promise.
    */
    const key = _queryCacheKey(field, query);
    if (!this._pendingLoad[field][key]) {
      this._pendingLoad[field][key] = fetchFn(field, query);
      try {
        await this._pendingLoad[field][key];
      } finally {
        delete this._pendingLoad[field][key];
      }
    }
    return this._pendingLoad[field][key];
  }

  // eslint-disable-next-line class-methods-use-this -- make sure subclass implements
  async _doLoad() {
    /* protect against bugs in subclass */
    throw new Error("subclass failed to implement _doLoad");
  }

  /**
   ** Garbage collection of annomatrix cache to manage memory use.
   **/

  /*
  These callbacks implement a GC policy for the cache.  Background:

    * For the Loader (base) annomatrix, re-filling the cache is expensive as
      it requires an HTTP fetch.
    * user-defined / writable columns must not be GC'ed as they may be
      still pending a save/commit.
    * For views, cost is less and (roughly) proportional with nObs
    * obs, var and emb do not grow without bounds, and are needed constantly 
      for rendering.  
      a) There is no upside to GC'ing these in the base (loader)
      b) The undo/redo cache can hold a large number in views, which is worht GC'ing
    * X is often much larger than memory, and the UI allows add/del from
      this.  Most of the GC potential is here in both the base and views.

  Current policy:
    * if in active use ("hot") do not GC obs, var or emb.
    * never, ever GC writable obs columns
    * For base/loader set a numeric limit on maximum X column count
    * For views, apply a fixed limit to the number of columns cached in any field.
      Limit will be lower if not hot.

  To be effective, the GC callback needs to be invoked from the undo/redo code,
  as much of the cache is pinned by that data structure.
  */
  _gcField(field, isHot, pinnedColumns) {
    const maxColumns = isHot ? 256 : 10; // maybe to aggessive?

    if (this[field].length < maxColumns) return;

    const gcInfo = this._gcInfo[field];
    const candidates = Array.from(gcInfo.keys()).filter(
      (col) => !pinnedColumns.includes(col)
    );

    const excessCount = candidates.length + pinnedColumns.length - maxColumns;
    if (excessCount > 0) {
      const toDrop = candidates.slice(0, excessCount);
      console.log(`dropping from ${field} hot:${isHot}, columns [${toDrop.join(', ')}]`);
      this[field] = toDrop.reduce((df, col) => df.dropCol(col), this[field]);
      toDrop.forEach((col) => gcInfo.delete(col));
    }
  }

  _gcFetchCleanup(field, pinnedColumns) {
    /*
    Called during data load/fetch.  By definition, this is 'hot', so we
    only want to gc X.
    */
    if (field === "X") {
      this._gcField(
        field,
        true,
        pinnedColumns.concat(_getWritableColumns(this.schema, field))
      );
    }
  }

  _gc(hints) {
    /*
    Called from middleware, or elsewhere.  isHot is true if we are in the active store, 
    or false if we are in some other context (eg, history state).
    */
    const { isHot } = hints;
    const candidateFields = isHot ? ["X"] : ["X", "emb", "var", "obs"];
    candidateFields.forEach((field) =>
      this._gcField(field, isHot, _getWritableColumns(this.schema, field))
    );
  }

  _gcUpdateStats(field, dataframe) {
    /*
    called each time a query is performed, allowing the gc to update any bookkeeping
    information.  Currently, this is just a simple last-fetched timestamp, stored
    in a Map.

    Map objects preserve order of insertion. This is leveraged as a cheap way to
    do LRU, by removing and re-inserting keys.  IMPORTANT: the cleanup code assumes
    the map insertion order is least-recently-used first.
    */
    const cols = dataframe.colIndex.labels();
    const gcInfo = this._gcInfo[field];
    const now = Date.now();
    cols.forEach((c) => {
      gcInfo.delete(c);
      gcInfo.set(c, now);
    });
  }
}

/*
private utility functions below
*/

function _queryCacheKey(field, query) {
  if (typeof query === "object") {
    const { field: queryField, column: queryColumn, value: queryValue } = query;
    return `${field}/${queryField}/${queryColumn}/${queryValue}`;
  }
  return `${field}/${query}`;
}