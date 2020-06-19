/*
Base matrix, which proxies for server.
*/
import clone from "./clone";
import { doBinaryRequest, dubEncURIComp } from "./fetchHelpers";
import { matrixFBSToDataframe } from "../stateManager/matrix";
import { getColumnSchema } from "./schema";
import {
  addObsAnnoColumn,
  removeObsAnnoColumn,
  addObsAnnoCategory,
  removeObsAnnoCategory,
} from "../stateManager/schemaHelpers";
import { isArrayOrTypedArray } from "../typeHelpers";
import { whereCacheCreate } from "./whereCache";
import AnnoMatrix from "./annoMatrix";

export default class AnnoMatrixLoader extends AnnoMatrix {
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
    /*
    Add a new category (aka label) to the schema for an obs column.
    */
    const colSchema = getColumnSchema(this.schema, "obs", col);
    if (!colSchema.writable) {
      throw new Error("Not a writable obs column");
    }

    const o = clone(this);
    o.schema = addObsAnnoCategory(this.schema, col, category);
    return o;
  }

  async removeObsAnnoCategory(col, category, unassignedCategory) {
    /*
    Remove a single "category" (aka "label") from the data & schema of an obs column.
    */
    const colSchema = getColumnSchema(this.schema, "obs", col);
    if (!colSchema.writable) {
      throw new Error("Not a writable obs column");
    }

    const o = await this.resetObsColumnValues(
      col,
      category,
      unassignedCategory
    );
    o.schema = removeObsAnnoCategory(this.schema, col, category);
    return o;
  }

  dropObsColumn(col) {
    /*
		drop column from field
		*/
    const colSchema = getColumnSchema(this.schema, "obs", col);
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
    const col = colSchema.name;
    if (getColumnSchema(this.schema, "obs", col) || this.obs.hasCol(col)) {
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
    const oldColSchema = getColumnSchema(this.schema, "obs", oldCol);
    if (!oldColSchema.writable) {
      throw new Error("Not a writable obs column");
    }

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
    const colSchema = getColumnSchema(this.schema, "obs", col);
    if (!colSchema.writable) {
      throw new Error("Not a writable obs column");
    }

    // ensure that we have the data in cache before we manipulate it
    await this.fetch("obs", col);
    if (!this.obs.hasCol(col))
      throw new Error("Internal error - user annotation data missing");

    const rowIndices = this.rowIndex.getOffsets(rowLabels);
    const data = this.obs.col(col).asArray().slice();
    for (let i = 0, len = rowIndices.length; i < len; i += 1) {
      data[rowIndices[i]] = value;
    }

    const o = clone(this);
    o.obs = this.obs.replaceColData(col, data);
    const { categories } = colSchema;
    if (categories && categories.indexOf(value) === -1) {
      o.schema = addObsAnnoCategory(this.schema, col, value);
    }
    return o;
  }

  async resetObsColumnValues(col, oldValue, newValue) {
    /*
    Set all rows with value 'oldValue' to 'newValue'
    */
    const colSchema = getColumnSchema(this.schema, "obs", col);
    if (!colSchema.writable) {
      throw new Error("Not a writable obs column");
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
    if (categories && categories.indexOf(newValue) === -1) {
      o.schema = addObsAnnoCategory(this.schema, col, newValue);
    }
    return o;
  }
}

/*
Utility functions below
*/

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
