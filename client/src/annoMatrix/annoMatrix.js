import {
  Dataframe,
  IdentityInt32Index,
  dataframeMemo,
} from "../util/dataframe";
import {
  _getColumnDimensionNames,
  _getColumnSchema,
  _schemaColumns,
  _getWritableColumns,
} from "./schema";
import { indexEntireSchema } from "../util/stateManager/schemaHelpers";
import { _whereCacheGet, _whereCacheMerge } from "./whereCache";
import _shallowClone from "./clone";
import { _queryValidate, _queryCacheKey } from "./query";

const _dataframeCache = dataframeMemo(128);

export default class AnnoMatrix {
  /*
  Abstract base class for all AnnoMatrix objects.  This class provides a proxy
  to the annotated matrix data authoritatively served by the server/back-end.

  AnnoMatrix instances are immutable, meaning that their schema and dimensionality
  will not change, and simple object equality can be used to detect structural
  changes. The actual data is cached, and not guaranteed to be present -- any
  request to access data must be resolved by a fetch() call, which is async, and
  may involve a server round-trip.

  Guarantees made by the immutabilty, ie, any of these can be detected by
  simple annoMatrix compare:
    * schema is the same, including all fields and columns
    * dimensionality is the same (nObs, nVar)
    * data mapping/transformation, such as clipping, are the same

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
    return ["obs", "var", "emb", "X", "layers"];
  }

  constructor(schema, nObs, nVar, rowIndex = null) {
    /*
    Private constructor - this is an abstract base class.  Do not use.
    */

    /*
    Public instance fields:
      * schema - the matrix schema.  IMPORTANT: always the entire schema, for the
        base (unfiltered, unclipped, unsubset) annotated matrix, as the server
        presents it.
      * nObs, nVar - size of each dimension.  These will accurately reflect the
        size of the current annoMatrix view.  For example, if you subset the view,
        the nObs will be smaller.
      * rowIndex - a rowIndex shared by all data on this view (ie, the list of cells).
        The row index labels are as defined by the base dataset from the server.
      * isView - true if this is a view, false if not.
      * viewOf - pointer to parent annomatrix if a view, undefined/null if not a view.
      * userFlags - container for any additional state a user of this API wants to hang
        off of an annoMatrix, and have propagated by the (shallow) cloning protocol.
    */
    this.schema = indexEntireSchema(schema);
    this.nObs = nObs;
    this.nVar = nVar;
    this.rowIndex = rowIndex || new IdentityInt32Index(nObs);
    this.isView = false;
    this.viewOf = undefined;
    this.userFlags = {};

    /*
		Private instance variables.

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
    /*
    Return array of column names in the field.  ONLY supported on the
    obs, var and emb fields.  X currently unimplemented and will throw.

    For exmaple:

      annoMatrix.getMatrixColumns("obs") -> ["louvain", "n_genes"]
    */
    return _schemaColumns(this.schema, field);
  }

  // eslint-disable-next-line class-methods-use-this -- need to be able to call this on instances
  getMatrixFields() {
    /*
    Return array of fields in this annoMatrix.  Currently hard-wired to
    return:  ["X", "obs", "var", "emb"].

    These are the fields from data may be requested.
    */
    return AnnoMatrix.fields();
  }

  getColumnSchema(field, col) {
    /*
    Return the schema for the field & column ,eg,

      anonMatrix.getColumnSchema("obs", "n_genes") -> { type: "int32", name: "n_genes" }

    This is identical to the information in the annoMatrix.schema
    instance variable.
    */
    return _getColumnSchema(this.schema, field, col);
  }

  getColumnDimensions(field, col) {
    /*
    Return the dimensions on this field / column.  For most fields, which are 1D,
    this just return the column name.  Multi-dimensional columns, such as embeddings,
    will return >1 name.

    Examples:

      getColumnDimensions("obs", "louvain") -> ["louvain"]
      getColumnDimensions("emb", "umap") -> ["umap_0", "umap_1"]

    */
    return _getColumnDimensionNames(this.schema, field, col);
  }

  /**
   ** General utility methods
   **/
  base() {
    /*
    return the base of view, or `this` if not a view.
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

		Returns a Promise for the query result, which will resolve to a dataframe.

		Field must be one of the matrix fields: 'obs', 'var', 'X', 'emb'.  Value
		represents the underlying object upon which the query is occurring.

		Query is one of:
			* a string, representing a single column name from the field, eg,
				"n_genes"
			* an object, containing an "value" query (see below).
			* an array, containing one or more of the above.

		Columns may have more than one dimension, and all will be fetched
		and returned together.  This is most commonly seen in an embedding,
    which usually has two dimensions.

		A value query allows for fetching based upon the value in another 
		field/column, similar to a join.  Currently only supported on the var
		dimension, allowing query of X columns by var value (eg, gene name)

		Examples:

    1. Fetch the "n_genes" column the "obs":

			const df = await fetch("obs", "n_genes")
      console.log("Largest number of genes is: ", df.summarize().max);

    2. Fetch two separate columns from obs.  Returns a single dataframe containing
       the columns:

			const df = await fetch("obs", ["n_genes", "louvain"])
      console.log("Cell 0 has category: ", df.at(0, "louvain"));

    3. Fetch an entire X (expression counts) column that has a var annotation
       value "TYMP" in the var index.

			fetch("X", { 
				where: {
          field: "var", column: this.schema.annotations.var.index, value: "TYMP"
        }
			})

      In AnnData & Pandas DataFrame API, this is equivalent to:
        adata.X[:, adata.var.index.get_loc("SUMO3")]

		The value query is a recodification and subset of the server REST API 
		value filter JSON.  Range queries and multiple filters are not currently
		supported.

		*/
    return this._fetch(field, q);
  }

  prefetch(field, q) {
    /*
		Start a data fetch & cache fill.  Identical to fetch() except it does
		not return a value.

    Primary use is to being a cache load as early as is possible, reducing
    overall component rendering latency.
		*/
    this._fetch(field, q);
    return undefined;
  }

  /**
   ** Save / mutate interfaces - manipulation of "writable" OBS annotations.
   **
   ** These are all present to support client-side creation of OBS annotations, aka
   ** "user annotations".
   **
   ** They implement common manipulations to the AnnoMatrix, maintaining the
   ** norma guarantees around correctness of public API, eg,
   **   - schema will be correct, including the "writable" attribute
   **   - fetch() will return the latest data, even from views
   **   - immutability guranteeds
   **
   ** As most of these interfaces mutate the annoMatrix, they return a new
   ** annoMatrix
   **
   ** The actual implementation is in the sub-classes, which MUST override these.
   **/

  // eslint-disable-next-line class-methods-use-this, no-unused-vars -- make sure subclass implements
  addObsAnnoCategory(col, category) {
    /*
    Add a new category value (aka "label") to a writable obs column, and return the new AnnoMatrix.
    Typical use is to add a new user-created label to a user-created obs categorical
    annotation.

    Will throw column does not exist or is not writable.

    Example:

      addObsAnnoCategory("my cell type", "left toenail") -> AnnoMatrix

    */
    _subclassResponsibility();
  }

  // eslint-disable-next-line class-methods-use-this, no-unused-vars -- make sure subclass implements
  async removeObsAnnoCategory(col, category, unassignedCategory) {
    /*
    Remove a category value from an obs column, reassign any obs having that value
    to the 'unassignedCategory' value, and return a promise for a new AnnoMatrix.
    Typical use is to remove a user-created label from a user-created obs categorical
    annotation.

    Will throw column does not exist or is not writable.

    An `unassignedCategory` value must be provided, for assignment to any obs/cells
    that had the now-delete category label as their value.

    Example:
      await removeObsAnnoCategory("my-tissue-type", "right earlobe", "unassigned") -> AnnoMatrix

    NOTE: method is async as it may need to fetch data to provide the reassignment.
    */
    _subclassResponsibility();
  }

  // eslint-disable-next-line class-methods-use-this, no-unused-vars -- make sure subclass implements
  dropObsColumn(col) {
    /*
    Drop an entire writable column, eg a user-created obs annotation.  Typical use
    is to provide the "Delete Category" implementation.  Returns the new AnnoMatrix.
    Will throw if not a writable annotation.

    Will throw column does not exist or is not writable.

    Example:

      dropObsColumn("old annotations") ->  AnnoMatrix
    */
    _subclassResponsibility();
  }

  // eslint-disable-next-line class-methods-use-this, no-unused-vars -- make sure subclass implements
  addObsColumn(colSchema, Ctor, value) {
    /*
    Add a new writable OBS annotation column, with the caller-specified schema, initial value
    type and value.

    Value may be any one of:
      * an array of values
      * a primitive type, including null or undefined.
    If an array, length must be the same as 'this.nObs', and constructor must equal 'Ctor'.
    If a primitive, 'Ctor' will be used to create the initial value, which will be filled
    with 'value'.

    Throws if the name specified in 'colSchema' duplicates an existing obs column.

    Returns a new AnnoMatrix.

    Examples:

      addObsColumn(
        { name: "foo", type: "categorical", categories: "unassigned" },
        Array,
        "unassigned"
      ) -> AnnoMatrix

    */
    _subclassResponsibility();
  }

  // eslint-disable-next-line class-methods-use-this, no-unused-vars -- make sure subclass implements
  renameObsColumn(oldCol, newCol) {
    /*
    Rename the obs column 'oldCol' to have name 'newCol' and returns new AnnoMatrix.

    Will throw column does not exist or is not writable, or if 'newCol' is not unique.

    Example:

      renameObsColumn('cell type', 'old cell type') -> AnnoMatrix.

    */
    _subclassResponsibility();
  }

  // eslint-disable-next-line class-methods-use-this, no-unused-vars -- make sure subclass implements
  async setObsColumnValues(col, obsLabels, value) {
    /*
    Set all obs with label in array 'obsLabels' to have 'value'.  Typical use would be
    to set a group of cells to have a label on a user-created categorical anntoation
    (eg set all selected cells to have a label).

    NOTE: async method, as it may need to fetch.

    Will throw column does not exist or is not writable.

    Example:
      await setObsColmnValues("flavor", [383, 400], "tasty") -> AnnoMtarix

    */
    _subclassResponsibility();
  }

  // eslint-disable-next-line class-methods-use-this, no-unused-vars -- make sure subclass implements
  async resetObsColumnValues(col, oldValue, newValue) {
    /*
    Set by value - all elements in the column with value 'oldValue' are set to 'newValue'.
    Async method - returns a promise for a new AnnoMatrix.

    Typical use would be to set all labels of one value to another.

    Will throw column does not exist or is not writable.

    Example:
      await resetObsColumnValues("my notes", "good", "not-good") -> AnnoMatrix

    */
    _subclassResponsibility();
  }

  // eslint-disable-next-line class-methods-use-this, no-unused-vars -- make sure subclass implements
  addEmbedding(colSchema) {
    /*
    Add a new obs embedding to the AnnoMatrix, with provided schema.
    Returns a new annomatrix.

    Typical use will be to add a re-embedding that the server has calculated.

    Will throw if the column schema is invalid (eg, duplicate name).
    */
    _subclassResponsibility();
  }

  getCacheKeys(field, query) {
    /*
    Return cache keys for columns associated with this query.  May return
    [unknown] if no keys are known (ie, nothing is or was cached).
    */
    return _whereCacheGet(this._whereCache, this.schema, field, query);
  }

  /**
   ** Private interfaces below.
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
    queries.forEach(_queryValidate);

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
    const response = _dataframeCache(
      this._cache[field].subset(null, requestedCacheKeys)
    );
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
    _subclassResponsibility();
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
      b) The undo/redo cache can hold a large number in views, which is worth GC'ing
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
    const maxColumns = isHot ? 256 : 10; // maybe to aggressive?

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
      // helpful debugging - please leave in place.
      // console.log(
      //   `GC: dropping from ${field} hot:${isHot}, columns [${toDrop.join(
      //     ", "
      //   )}]`
      // );
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
    */
    const cols = dataframe.colIndex.labels();
    const { _gcInfo } = this;
    const now = Date.now();
    cols.forEach((c) => {
      _gcInfo.set(_columnCacheKey(field, c), now);
    });
  }

  /**
  Cloning subclass protocol - we rely in cloning to preserve immutable
  semantics while not causing races or other side effects in internal
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
function _columnCacheKey(field, column) {
  return `${field}/${column}`;
}

function _subclassResponsibility() {
  /* protect against bugs in subclass */
  throw new Error("subclass failed to implement required method");
}
