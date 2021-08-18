import { doBinaryRequest, doFetch } from "./fetchHelpers";
import { matrixFBSToDataframe } from "../util/stateManager/matrix";
import { _getColumnSchema } from "./schema";
import {
  addObsAnnoCategory,
  addObsAnnoColumn,
  addObsLayout,
  removeObsAnnoCategory,
  removeObsAnnoColumn,
} from "../util/stateManager/schemaHelpers";
import { isAnyArray } from "../common/types/arraytypes";
import { _whereCacheCreate, WhereCache } from "./whereCache";
import AnnoMatrix from "./annoMatrix";
import PromiseLimit from "../util/promiseLimit";
import {
  _expectComplexQuery,
  _expectSimpleQuery,
  _hashStringValues,
  _urlEncodeComplexQuery,
  _urlEncodeLabelQuery,
  ComplexQuery,
  Query,
} from "./query";
import {
  normalizeResponse,
  normalizeWritableCategoricalSchema,
} from "./normalize";
import {
  AnnotationColumnSchema,
  Field,
  EmbeddingSchema,
  RawSchema,
} from "../common/types/schema";
import {
  Dataframe,
  DataframeValue,
  DataframeValueArray,
  LabelType,
} from "../util/dataframe";

const promiseThrottle = new PromiseLimit(5);

export default class AnnoMatrixLoader extends AnnoMatrix {
  baseURL: string;

  /*
  AnnoMatrix implementation which proxies to HTTP server using the CXG REST API.
  Used as the base (non-view) instance.

  Public API is same as AnnoMatrix class (refer there for API description),
  with the addition of the constructor which bootstraps:

    new AnnoMatrixLoader(serverBaseURL, schema) -> instance

  */
  constructor(baseURL: string, schema: RawSchema) {
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
  addObsAnnoCategory(col: LabelType, category: string): AnnoMatrix {
    /*
    Add a new category (aka label) to the schema for an obs column.
    */
    const colSchema = _getColumnSchema(
      this.schema,
      Field.obs,
      col
    ) as AnnotationColumnSchema;
    _writableObsCategoryTypeCheck(colSchema); // throws on error

    const newAnnoMatrix = this._clone();
    newAnnoMatrix.schema = addObsAnnoCategory(this.schema, col, category);
    return newAnnoMatrix;
  }

  async removeObsAnnoCategory(
    col: LabelType,
    category: string,
    unassignedCategory: string
  ): Promise<AnnoMatrix> {
    /*
    Remove a single "category" (aka "label") from the data & schema of an obs column.
    */
    const colSchema = _getColumnSchema(
      this.schema,
      Field.obs,
      col
    ) as AnnotationColumnSchema;
    _writableObsCategoryTypeCheck(colSchema); // throws on error

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

  dropObsColumn(col: LabelType): AnnoMatrix {
    /*
		drop column from field
		*/
    const colSchema = _getColumnSchema(
      this.schema,
      Field.obs,
      col
    ) as AnnotationColumnSchema;
    _writableObsCheck(colSchema); // throws on error

    const newAnnoMatrix = this._clone();
    newAnnoMatrix._cache.obs = this._cache.obs.dropCol(col);
    newAnnoMatrix.schema = removeObsAnnoColumn(this.schema, col);
    return newAnnoMatrix;
  }

  addObsColumn<T extends DataframeValueArray>(
    colSchema: AnnotationColumnSchema,
    Ctor: new (n: number) => T,
    value: T
  ): AnnoMatrix {
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
      _getColumnSchema(this.schema, Field.obs, colName) ||
      this._cache.obs.hasCol(colName)
    ) {
      throw new Error("column already exists");
    }

    const newAnnoMatrix = this._clone();
    let data;
    if (isAnyArray(value)) {
      if (value.constructor !== Ctor)
        throw new Error("Mismatched value array type");
      if (value.length !== this.nObs)
        throw new Error("Value array has incorrect length");
      data = value.slice();
    } else {
      data = new Ctor(this.nObs).fill(value);
    }
    newAnnoMatrix._cache.obs = this._cache.obs.withCol(colName, data);
    normalizeWritableCategoricalSchema(
      colSchema,
      newAnnoMatrix._cache.obs.col(colName)
    );
    newAnnoMatrix.schema = addObsAnnoColumn(this.schema, colName, colSchema);
    return newAnnoMatrix;
  }

  renameObsColumn(oldCol: LabelType, newCol: LabelType): AnnoMatrix {
    /*
    Rename the obs oldColName to newColName.  oldCol must be writable.
    */
    const oldColSchema = _getColumnSchema(
      this.schema,
      Field.obs,
      oldCol
    ) as AnnotationColumnSchema;
    _writableObsCheck(oldColSchema);
    const value = this._cache.obs.hasCol(oldCol)
      ? this._cache.obs.col(oldCol).asArray()
      : undefined;
    return this.dropObsColumn(oldCol).addObsColumn(
      {
        ...oldColSchema,
        // @ts-expect-error ts-migrate --- TODO revisit:
        // `name`: Type 'LabelType' is not assignable to type 'string'. Type 'number' is not assignable to type 'string'.
        name: newCol,
      },
      // @ts-expect-error ts-migrate --- TODO revisit:
      // `value`: Object is possibly 'undefined'.
      value.constructor,
      value
    );
  }

  async setObsColumnValues(
    col: LabelType,
    rowLabels: Int32Array,
    value: DataframeValue
  ): Promise<AnnoMatrix> {
    /*
		Set all rows identified by rowLabels to value.
		*/
    const colSchema = _getColumnSchema(
      this.schema,
      Field.obs,
      col
    ) as AnnotationColumnSchema;
    _writableObsCategoryTypeCheck(colSchema); // throws on error

    // ensure that we have the data in cache before we manipulate it
    await this.fetch(Field.obs, col);
    if (!this._cache.obs.hasCol(col))
      throw new Error("Internal error - user annotation data missing");

    const rowIndices = this.rowIndex.getOffsets(rowLabels);
    const data = this._cache.obs.col(col).asArray().slice();
    for (let i = 0, len = rowIndices.length; i < len; i += 1) {
      const idx = rowIndices[i];
      if (idx === -1) throw new Error("Unknown row label");
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

  async resetObsColumnValues<T extends DataframeValue>(
    col: LabelType,
    oldValue: T,
    newValue: T
  ): Promise<AnnoMatrix> {
    /*
    Set all rows with value 'oldValue' to 'newValue'.
    */
    const colSchema = _getColumnSchema(
      this.schema,
      Field.obs,
      col
    ) as AnnotationColumnSchema;
    _writableObsCategoryTypeCheck(colSchema); // throws on error

    // @ts-expect-error ts-migrate --- TODO revisit:
    // `colSchema.categories`: Object is possibly 'undefined'.
    if (!colSchema.categories.includes(oldValue)) {
      throw new Error("unknown category");
    }

    // ensure that we have the data in cache before we manipulate it
    await this.fetch(Field.obs, col);
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

  addEmbedding(colSchema: EmbeddingSchema): AnnoMatrix {
    /*
    add new layout to the obs embeddings
    */
    const { name: colName } = colSchema;
    if (_getColumnSchema(this.schema, Field.emb, colName)) {
      throw new Error("column already exists");
    }

    const newAnnoMatrix = this._clone();
    newAnnoMatrix.schema = addObsLayout(this.schema, colSchema);
    return newAnnoMatrix;
  }

  /**
   ** Private below
   **/
  async _doLoad(
    field: Field,
    query: Query
  ): Promise<[WhereCache | null, Dataframe]> {
    /*
    _doLoad - evaluates the query against the field. Returns:
      * whereCache update: column query map mapping the query to the column labels
      * Dataframe containing the new columns (one per dimension)
    */
    let doRequest;
    let priority = 10; // default fetch priority

    switch (field) {
      case "obs":
      case "var": {
        doRequest = _obsOrVarLoader(this.baseURL, field, query);
        break;
      }
      case "X": {
        doRequest = _XLoader(this.baseURL, field, query);
        break;
      }
      case "emb": {
        doRequest = _embLoader(this.baseURL, field, query);
        priority = 0; // high prio load for embeddings
        break;
      }
      default:
        throw new Error("Unknown field name");
    }
    const buffer = await promiseThrottle.priorityAdd(priority, doRequest);
    // @ts-expect-error --- TODO revisit:
    // `buffer`:  Argument of type 'unknown' is not assignable to parameter of type 'ArrayBuffer | ArrayBuffer[]'. Type 'unknown' is not assignable to type 'ArrayBuffer[]'.
    let result = matrixFBSToDataframe(buffer);
    if (!result || result.isEmpty()) throw Error("Unknown field/col");

    const whereCacheUpdate = _whereCacheCreate(
      field,
      query,
      result.colIndex.labels()
    );

    result = normalizeResponse(field, this.schema, result);

    return [whereCacheUpdate, result];
  }
}

/*
Utility functions below
*/

function _writableObsCheck(obsColSchema: AnnotationColumnSchema): void {
  if (!obsColSchema?.writable) {
    throw new Error("Unknown or readonly obs column");
  }
}

function _writableObsCategoryTypeCheck(
  obsColSchema: AnnotationColumnSchema
): void {
  _writableObsCheck(obsColSchema);
  if (obsColSchema.type !== "categorical") {
    throw new Error("column must be categorical");
  }
}

function _embLoader(
  baseURL: string,
  _field: Field,
  query: Query
): () => Promise<ArrayBuffer> {
  _expectSimpleQuery(query);

  const urlBase = `${baseURL}layout/obs`;
  const urlQuery = _urlEncodeLabelQuery("layout-name", query);
  const url = `${urlBase}?${urlQuery}`;
  return () => doBinaryRequest(url);
}

function _obsOrVarLoader(
  baseURL: string,
  field: Field,
  query: Query
): () => Promise<ArrayBuffer> {
  _expectSimpleQuery(query);

  const urlBase = `${baseURL}annotations/${field}`;
  const urlQuery = _urlEncodeLabelQuery("annotation-name", query);
  const url = `${urlBase}?${urlQuery}`;
  return () => doBinaryRequest(url);
}

function _XLoader(
  baseURL: string,
  _field: Field,
  query: Query
): () => Promise<ArrayBuffer> {
  _expectComplexQuery(query);

  // Casting here as query is validated to be complex in _expectComplexQuery above.
  const complexQuery = query as ComplexQuery;

  if ("where" in complexQuery) {
    const urlBase = `${baseURL}data/var`;
    const urlQuery = _urlEncodeComplexQuery(complexQuery);
    const url = `${urlBase}?${urlQuery}`;
    return () => doBinaryRequest(url);
  }

  if ("summarize" in complexQuery) {
    const urlBase = `${baseURL}summarize/var`;
    const urlQuery = _urlEncodeComplexQuery(complexQuery);

    if (urlBase.length + urlQuery.length < 2000) {
      const url = `${urlBase}?${urlQuery}`;
      return () => doBinaryRequest(url);
    }

    const url = `${urlBase}?key=${_hashStringValues([urlQuery])}`;
    return async () => {
      const res = await doFetch(url, {
        method: "POST",
        body: urlQuery,
        headers: new Headers({
          Accept: "application/octet-stream",
          "Content-Type": "application/x-www-form-urlencoded",
        }),
      });
      return res.arrayBuffer();
    };
  }

  throw new Error("Unknown query structure");
}
