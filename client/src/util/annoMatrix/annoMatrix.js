import { Dataframe, IdentityInt32Index } from "../dataframe";
import {
  getColumnDimensionNames,
  getColumnSchema,
  schemaColumns,
} from "./schema";
import { fetchResult } from "./fetchHelpers";
import { indexEntireSchema } from "../stateManager/schemaHelpers";
import { whereCacheGet, whereCacheMerge } from "./whereCache";

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
Utility functions below
*/

function queryCacheKey(field, query) {
  if (typeof query === "object") {
    const { field: queryField, column: queryColumn, value: queryValue } = query;
    return `${field}/${queryField}/${queryColumn}/${queryValue}`;
  }
  return `${field}/${query}`;
}
