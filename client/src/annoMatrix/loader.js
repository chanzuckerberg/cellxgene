/*
Base matrix, which proxies for server.
*/
import clone from "./clone";
import { doBinaryRequest, _dubEncURIComp } from "./fetchHelpers";
import { matrixFBSToDataframe } from "../util/stateManager/matrix";
import { _getColumnSchema, _normalizeCategoricalSchema } from "./schema";
import {
  addObsAnnoColumn,
  removeObsAnnoColumn,
  addObsAnnoCategory,
  removeObsAnnoCategory,
} from "../util/stateManager/schemaHelpers";
import { isArrayOrTypedArray } from "../util/typeHelpers";
import { _whereCacheCreate } from "./whereCache";
import AnnoMatrix from "./annoMatrix";

export default class AnnoMatrixLoader extends AnnoMatrix {
  constructor(baseURL, schema) {
    const { nObs, nVar } = schema.dataframe;
    super(schema, nObs, nVar);

    if (baseURL[baseURL.length - 1] !== "/") {
      // must have trailing slash
      baseURL += "/";
    }
    this.baseURL = baseURL;
  }

  /**
   ** Public
   **/
  addObsAnnoCategory(col, category) {
    /*
    Add a new category (aka label) to the schema for an obs column.
    */
    const colSchema = _getColumnSchema(this.schema, "obs", col);
    _writableCategoryTypeCheck(colSchema); // throws on error

    const o = clone(this);
    o.schema = addObsAnnoCategory(this.schema, col, category);
    return o;
  }

  async removeObsAnnoCategory(col, category, unassignedCategory) {
    /*
    Remove a single "category" (aka "label") from the data & schema of an obs column.
    */
    const colSchema = _getColumnSchema(this.schema, "obs", col);
    _writableCategoryTypeCheck(colSchema); // throws on error

    const o = await this.resetObsColumnValues(
      col,
      category,
      unassignedCategory
    );
    o.schema = removeObsAnnoCategory(o.schema, col, category);
    return o;
  }

  dropObsColumn(col) {
    /*
		drop column from field
		*/
    const colSchema = _getColumnSchema(this.schema, "obs", col);
    _writableCheck(colSchema); // throws on error

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
    colSchema.writable = true;
    const col = colSchema.name;
    if (_getColumnSchema(this.schema, "obs", col) || this.obs.hasCol(col)) {
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
    const oldColSchema = _getColumnSchema(this.schema, "obs", oldCol);
    _writableCheck(oldColSchema); // throws on error

    const value = this.obs.hasCol(oldCol)
      ? this.obs.col(oldCol).asArray()
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
    if (!this.obs.hasCol(col))
      throw new Error("Internal error - user annotation data missing");

    const rowIndices = this.rowIndex.getOffsets(rowLabels);
    const data = this.obs.col(col).asArray().slice();
    for (let i = 0, len = rowIndices.length; i < len; i += 1) {
      const idx = rowIndices[i];
      if (idx === undefined) throw new Error("Unknown row label");
      data[idx] = value;
    }

    const o = clone(this);
    o.obs = this.obs.replaceColData(col, data);
    const { categories } = colSchema;
    if (!categories?.includes(value)) {
      o.schema = addObsAnnoCategory(this.schema, col, value);
    }
    return o;
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
    if (!this.obs.hasCol(col))
      throw new Error("Internal error - user annotation data missing");

    const data = this.obs.col(col).asArray().slice();
    for (let i = 0, l = data.length; i < l; i += 1) {
      if (data[i] === oldValue) data[i] = newValue;
    }

    const o = clone(this);
    o.obs = this.obs.replaceColData(col, data);
    const { categories } = colSchema;
    if (!categories?.includes(newValue)) {
      o.schema = addObsAnnoCategory(this.schema, col, newValue);
    }
    return o;
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
        break;
      }
      default:
        throw new Error("Unknown field name");
    }

    const url = `${urlBase}?${urlQuery}`;
    const buffer = await doBinaryRequest(url);
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
