import { doBinaryRequest, _dubEncURIComp } from "./fetchHelpers";
import { matrixFBSToDataframe } from "../util/stateManager/matrix";
import { _getColumnSchema, _normalizeCategoricalSchema } from "./schema";
import {
  addObsAnnoColumn,
  removeObsAnnoColumn,
  addObsAnnoCategory,
  removeObsAnnoCategory,
  addObsLayout,
} from "../util/stateManager/schemaHelpers";
import { isArrayOrTypedArray } from "../util/typeHelpers";
import { _whereCacheCreate } from "./whereCache";
import AnnoMatrix from "./annoMatrix";
import PromiseLimit from "../util/promiseLimit";

const promiseThrottle = new PromiseLimit(5);

export default class AnnoMatrixLoader extends AnnoMatrix {
  /*
  AnnoMatrix implementation which proxies to HTTP server using the CXG REST API.
  Used as the base (non-view) instance.

  Public API is same as AnnoMatrix class (refer there for API description),
  with the addition of the constructor which bootstraps:

    new AnnoMatrixLoader(serverBaseURL, schema) -> instance

  */
  constructor(baseURL, schema) {
    const { nObs, nVar } = schema.dataframe;
    super(schema, nObs, nVar);

    if (baseURL[baseURL.length - 1] !== "/") {
      // must have trailing slash
      baseURL += "/";
    }
    this.baseURL = baseURL;
    Object.seal(this);
  }

  /**
   ** Public.  API described in base class.
   **/
  addObsAnnoCategory(col, category) {
    /*
    Add a new category (aka label) to the schema for an obs column.
    */
    const colSchema = _getColumnSchema(this.schema, "obs", col);
    _writableCategoryTypeCheck(colSchema); // throws on error

    const newAnnoMatrix = this._clone();
    newAnnoMatrix.schema = addObsAnnoCategory(this.schema, col, category);
    return newAnnoMatrix;
  }

  async removeObsAnnoCategory(col, category, unassignedCategory) {
    /*
    Remove a single "category" (aka "label") from the data & schema of an obs column.
    */
    const colSchema = _getColumnSchema(this.schema, "obs", col);
    _writableCategoryTypeCheck(colSchema); // throws on error

    const newAnnoMatrix = await this.resetObsColumnValues(
      col,
      category,
      unassignedCategory
    );
    newAnnoMatrix.schema = removeObsAnnoCategory(
      newAnnoMatrix.schema,
      col,
      category
    );
    return newAnnoMatrix;
  }

  dropObsColumn(col) {
    /*
		drop column from field
		*/
    const colSchema = _getColumnSchema(this.schema, "obs", col);
    _writableCheck(colSchema); // throws on error

    const newAnnoMatrix = this._clone();
    newAnnoMatrix._cache.obs = this._cache.obs.dropCol(col);
    newAnnoMatrix.schema = removeObsAnnoColumn(this.schema, col);
    return newAnnoMatrix;
  }

  addObsColumn(colSchema, Ctor, value) {
    /*
		add a column to field, initializing with value.  Value may 
    be one of:
      * an array of values
      * a primitive type, including null or undefined.
    If an array, it must be of same size as nObs and same type as Ctor
		*/
    colSchema.writable = true;
    const colName = colSchema.name;
    if (
      _getColumnSchema(this.schema, "obs", colName) ||
      this._cache.obs.hasCol(colName)
    ) {
      throw new Error("column already exists");
    }

    const newAnnoMatrix = this._clone();
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
    newAnnoMatrix._cache.obs = this._cache.obs.withCol(colName, data);
    _normalizeCategoricalSchema(
      colSchema,
      newAnnoMatrix._cache.obs.col(colName)
    );
    newAnnoMatrix.schema = addObsAnnoColumn(this.schema, colName, colSchema);
    return newAnnoMatrix;
  }

  renameObsColumn(oldCol, newCol) {
    /*
    Rename the obs oldColName to newColName.  oldCol must be writable.
    */
    const oldColSchema = _getColumnSchema(this.schema, "obs", oldCol);
    _writableCheck(oldColSchema); // throws on error

    const value = this._cache.obs.hasCol(oldCol)
      ? this._cache.obs.col(oldCol).asArray()
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

  async setObsColumnValues(col, rowLabels, value) {
    /*
		Set all rows identified by rowLabels to value.
		*/
    const colSchema = _getColumnSchema(this.schema, "obs", col);
    _writableCategoryTypeCheck(colSchema); // throws on error

    // ensure that we have the data in cache before we manipulate it
    await this.fetch("obs", col);
    if (!this._cache.obs.hasCol(col))
      throw new Error("Internal error - user annotation data missing");

    const rowIndices = this.rowIndex.getOffsets(rowLabels);
    const data = this._cache.obs.col(col).asArray().slice();
    for (let i = 0, len = rowIndices.length; i < len; i += 1) {
      const idx = rowIndices[i];
      if (idx === undefined) throw new Error("Unknown row label");
      data[idx] = value;
    }

    const newAnnoMatrix = this._clone();
    newAnnoMatrix._cache.obs = this._cache.obs.replaceColData(col, data);
    const { categories } = colSchema;
    if (!categories?.includes(value)) {
      newAnnoMatrix.schema = addObsAnnoCategory(this.schema, col, value);
    }
    return newAnnoMatrix;
  }

  async resetObsColumnValues(col, oldValue, newValue) {
    /*
    Set all rows with value 'oldValue' to 'newValue'.
    */
    const colSchema = _getColumnSchema(this.schema, "obs", col);
    _writableCategoryTypeCheck(colSchema); // throws on error

    if (!colSchema.categories.includes(oldValue)) {
      throw new Error("unknown category");
    }

    // ensure that we have the data in cache before we manipulate it
    await this.fetch("obs", col);
    if (!this._cache.obs.hasCol(col))
      throw new Error("Internal error - user annotation data missing");

    const data = this._cache.obs.col(col).asArray().slice();
    for (let i = 0, l = data.length; i < l; i += 1) {
      if (data[i] === oldValue) data[i] = newValue;
    }

    const newAnnoMatrix = this._clone();
    newAnnoMatrix._cache.obs = this._cache.obs.replaceColData(col, data);
    const { categories } = colSchema;
    if (!categories?.includes(newValue)) {
      newAnnoMatrix.schema = addObsAnnoCategory(this.schema, col, newValue);
    }
    return newAnnoMatrix;
  }

  addEmbedding(colSchema) {
    /*
    add new layout to the obs embeddings
    */
    const { name: colName } = colSchema;
    if (_getColumnSchema(this.schema, "emb", colName)) {
      throw new Error("column already exists");
    }

    const newAnnoMatrix = this._clone();
    newAnnoMatrix.schema = addObsLayout(this.schema, colSchema);
    return newAnnoMatrix;
  }

  /**
   ** Private below
   **/
  async _doLoad(field, query) {
    /*
    _doLoad - evaluates the query against the field. Returns:
      * whereCache update: column query map mapping the query to the column labels
      * Dataframe containing the new colums (one per dimension)
    */
    let urlQuery;
    let urlBase;
    let priority = 10; // default fetch priority

    switch (field) {
      case "obs":
      case "var": {
        urlBase = `${this.baseURL}annotations/${field}`;
        urlQuery = _encodeQuery("annotation-name", query);
        break;
      }
      case "X": {
        urlBase = `${this.baseURL}data/var`;
        urlQuery = _encodeQuery(undefined, query);
        break;
      }
      case "emb": {
        urlBase = `${this.baseURL}layout/obs`;
        urlQuery = _encodeQuery("layout-name", query);
        priority = 0; // high prio load for embeddings
        break;
      }
      default:
        throw new Error("Unknown field name");
    }

    const url = `${urlBase}?${urlQuery}`;
    const buffer = await promiseThrottle.priorityAdd(
      priority,
      doBinaryRequest,
      url
    );
    const result = matrixFBSToDataframe(buffer);
    if (!result || result.isEmpty()) throw Error("Unknown field/col");

    const whereCacheUpdate = _whereCacheCreate(
      field,
      query,
      result.colIndex.labels()
    );

    if (field === "obs") {
      /* cough, cough - see comment on method */
      _normalizeCategoricalSchema(
        this.schema.annotations.obsByName[query],
        result.col(query)
      );
    }

    return [whereCacheUpdate, result];
  }
}

/*
Utility functions below
*/

function _encodeQuery(colKey, q) {
  if (typeof q === "object") {
    const { field: queryField, column: queryColumn, value: queryValue } = q;
    return `${_dubEncURIComp(queryField)}:${_dubEncURIComp(
      queryColumn
    )}=${_dubEncURIComp(queryValue)}`;
  }
  if (!colKey) throw new Error("Unsupported query by name");
  return `${colKey}=${encodeURIComponent(q)}`;
}

function _writableCheck(colSchema) {
  if (!colSchema?.writable) {
    throw new Error("Unknown or readonly obs column");
  }
}

function _writableCategoryTypeCheck(colSchema) {
  _writableCheck(colSchema);
  if (colSchema.type !== "categorical") {
    throw new Error("column must be categorical");
  }
}
