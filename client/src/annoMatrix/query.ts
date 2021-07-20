import sha1 from "sha1";
import { _dubEncURIComp } from "./fetchHelpers";

/**
 * Query utilities, mostly for debugging support and validation.
 */

/**
 * Normalize & error check the query.
 * @param {object | string} query - the query
 * @returns {object | string} - the normalized query
 */
export function _queryValidate(query) {
  if (typeof query !== "object") return query;

  if (query.where && query.summarize)
    throw new Error("query may not specify both where and summarize");
  if (query.where) {
    const {
      field: queryField,
      column: queryColumn,
      value: queryValue,
    } = query.where;
    if (!queryField || !queryColumn || !queryValue)
      throw new Error("Incomplete where query");
    return query;
  }
  if (query.summarize) {
    const {
      field: queryField,
      column: queryColumn,
      values: queryValues,
    } = query.summarize;
    if (!queryField || !queryColumn || !queryValues)
      throw new Error("Incomplete where query");
    if (!Array.isArray(queryValues))
      throw new Error("Summarize query values must be an array");
    return query;
  }
  throw new Error("query must specify one of where or summarize");
}

export function _expectSimpleQuery(query) {
  if (typeof query === "object") throw new Error("expected simple query");
}

export function _expectComplexQuery(query) {
  if (typeof query !== "object") throw new Error("expected complex query");
}

/**
 * Generate a unique key which can be used to reference this query.
 *
 * @param {string} field
 * @param {string|object} query
 * @returns the key
 */
export function _queryCacheKey(field, query) {
  if (typeof query === "object") {
    // complex query
    if (query.where) {
      const {
        field: queryField,
        column: queryColumn,
        value: queryValue,
      } = query.where;
      return `${field}/${queryField}/${queryColumn}/${queryValue}`;
    }
    if (query.summarize) {
      const {
        method,
        field: queryField,
        column: queryColumn,
        values: queryValues,
      } = query.summarize;
      return `${field}/${method}/${queryField}/${queryColumn}/${queryValues.join(
        ","
      )}`;
    }
    throw new Error("Unrecognized complex query type");
  }

  // simple query
  return `${field}/${query}`;
}

function _urlEncodeWhereQuery(q) {
  const { field: queryField, column: queryColumn, value: queryValue } = q;
  return `${_dubEncURIComp(queryField)}:${_dubEncURIComp(
    queryColumn
  )}=${_dubEncURIComp(queryValue)}`;
}

function _urlEncodeSummarizeQuery(q) {
  const { method, field, column, values } = q;
  const filter = values
    .map((value) => _urlEncodeWhereQuery({ field, column, value }))
    .join("&");
  return `method=${method}&${filter}`;
}

export function _urlEncodeComplexQuery(q) {
  if (typeof q === "object") {
    if (q.where) {
      return _urlEncodeWhereQuery(q.where);
    }
    if (q.summarize) {
      return _urlEncodeSummarizeQuery(q.summarize);
    }
  }
  throw new Error("Unrecognized complex query type");
}

export function _urlEncodeLabelQuery(colKey, q) {
  if (!colKey) throw new Error("Unsupported query by name");
  if (typeof q !== "string") throw new Error("Query must be a simple label.");
  return `${colKey}=${encodeURIComponent(q)}`;
}

/**
 * Generate the column key the server will send us for this query.
 */
export function _hashStringValues(arrayOfString) {
  const hash = sha1(arrayOfString.join(""));
  return hash;
}
