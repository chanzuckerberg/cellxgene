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
import _shallowClone from "./clone";

export default class AnnoMatrix {
  /*
  Abstract base class for all AnnoMatrix objects.  This class provides a proxy
  to the annotated matrix data authoritatively served by the server/back-end.

  AnnoMatrix instances are immutable, meaning that their schema and dimensionality
  will not change, and simple object equality can be used to detect structural
  changes. The actual data is cached, and not guaranteed to be present -- any
  request to access data must be resolved by a fetch() call, which is async, and
  may involve a server round-trip.

  AnnoMatrixes also "stack" like filters, allowing for the construction of
  views which transform the data in some manner.

  The bootstrap class is AnnoMatrixLoader, which is the caching server proxy, and
  is bootstrapped with a API URL:
    new AnnoMatirx(url, schema) -> annoMatrix

  There are various "views", such as AnnoMatrixRowSubsetView, which provide
  the same interface but with a transformed view of the server data.  Utilities in
  viewCreators.js can be used to create these views:
      clip(annoMatrix, min, max) -> annoMatrix
      subset(annoMatrix, rowLabels) -> annoMatrix
  etc.
  */
  static fields() {
    /*
    return the fields present in the AnnoMatrix instance.
    */
    return ["obs", "var", "emb", "X"];
  }

  constructor(schema, nObs, nVar, rowIndex = null) {
    /*
    Private constructor - this is an abstract base class.
    */
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
    this._cache = {
      obs: Dataframe.empty(this.rowIndex),
      var: Dataframe.empty(this.rowIndex),
      emb: Dataframe.empty(this.rowIndex),
      X: Dataframe.empty(this.rowIndex),
    };
    this._pendingLoad = {
      obs: {},
      var: {},
      emb: {},
      X: {},
    };
    this._whereCache = {};
    this._gcInfo = new Map();
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
          (cacheKey) =>
            cacheKey !== undefined && this._cache[field].hasCol(cacheKey)
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
        (cacheKey) =>
          cacheKey === undefined || !this._cache[field].hasCol(cacheKey)
      )
    );

    /* load uncached queries */
    if (uncachedQueries.length > 0) {
      await Promise.all(
        uncachedQueries.map((query) =>
          this._getPendingLoad(field, query, async (_field, _query) => {
            /* fetch, then index.  _doLoad is subclass interface */
            const [whereCacheUpdate, df] = await this._doLoad(_field, _query);
            this._cache[_field] = this._cache[_field].withColsFrom(df);
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
    const response = this._cache[field].subset(null, requestedCacheKeys);
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

    const cache = this._cache[field];
    if (cache.colIndex.size() < maxColumns) return; // trivial rejection

    const candidates = cache.colIndex
      .labels()
      .filter((col) => !pinnedColumns.includes(col));

    const excessCount = candidates.length + pinnedColumns.length - maxColumns;
    if (excessCount > 0) {
      const { _gcInfo } = this;
      candidates.sort((a, b) => {
        let atime = _gcInfo.get(_columnCacheKey(field, a));
        if (atime === undefined) atime = 0;

        let btime = _gcInfo.get(_columnCacheKey(field, b));
        if (btime === undefined) btime = 0;

        return atime - btime;
      });

      const toDrop = candidates.slice(0, excessCount);
      console.log(
        `GC: dropping from ${field} hot:${isHot}, columns [${toDrop.join(
          ", "
        )}]`
      );
      this._cache[field] = toDrop.reduce(
        (df, col) => df.dropCol(col),
        this._cache[field]
      );
      toDrop.forEach((col) => _gcInfo.delete(_columnCacheKey(field, col)));
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
    const { _gcInfo } = this;
    const now = Date.now();
    cols.forEach((c) => {
      // gcInfo.delete(c);
      _gcInfo.set(_columnCacheKey(field, c), now);
    });
  }

  /**
  Cloning sublcass protocol - we rely in cloning to preserve immutable
  symantics while not causing races or other side effects in internal
  cache management.

  Subclasses must override _cloneDeeper() if they have state which requires 
  something other than a shallow copy.  Overrides MUST call super()._cloneDeepr(), 
  and return its result (after any required modification).  _cloneDeeper()
  will be called on the OLD object, with the NEW object as an argument.

  Do not override _clone();
  **/
  _cloneDeeper(clone) {
    clone._cache = _shallowClone(this._cache);
    clone._gcInfo = new Map();
    clone._pendingLoad = {
      obs: {},
      var: {},
      emb: {},
      X: {},
    };
    return clone;
  }

  _clone() {
    const clone = _shallowClone(this);
    this._cloneDeeper(clone);
    Object.seal(clone);
    return clone;
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

function _columnCacheKey(field, column) {
  return `${field}/${column}`;
}
