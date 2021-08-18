import sha1 from "sha1";
import { _dubEncURIComp } from "./fetchHelpers";
import { Field } from "../common/types/schema";
import { LabelType } from "../util/dataframe";

/**
 * Query utilities, mostly for debugging support and validation.
 */

export type ComplexQuery = SummarizeQuery | WhereQuery;

export type Query = LabelType | ComplexQuery;

interface SummarizeQuery {
  summarize: SummarizeQueryTerm;
}

interface SummarizeQueryTerm {
  column: string;
  field: string;
  method: string;
  values: string[];
}

interface WhereQuery {
  where: WhereQueryTerm;
}

interface WhereQueryTerm {
  column: string;
  field: string;
  value: string;
}

export function _expectSimpleQuery(query: Query): void {
  if (typeof query === "object") throw new Error("expected simple query");
}

/**
 * Normalize & error check the query.
 * @param {Query} query - the query
 * @returns {Query} - the normalized query
 */
export function _queryValidate(query: Query): Query {
  if (typeof query !== "object") return query;

  if ("where" in query && "summarize" in query)
    throw new Error("query may not specify both where and summarize");
  if ("where" in query) {
    const {
      field: queryField,
      column: queryColumn,
      value: queryValue,
    } = query.where;
    if (!queryField || !queryColumn || !queryValue)
      throw new Error("Incomplete where query");
    return query;
  }
  if ("summarize" in query) {
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

export function _expectComplexQuery(query: Query): void {
  if (typeof query !== "object") throw new Error("expected complex query");
}

/**
 * Generate a unique key which can be used to reference this query.
 *
 * @param {string} field
 * @param {string|object} query
 * @returns {string} the key
 */
export function _queryCacheKey(field: Field, query: Query): string {
  if (typeof query === "object") {
    // complex query
    if ("where" in query) {
      const {
        field: queryField,
        column: queryColumn,
        value: queryValue,
      } = query.where;
      return `${field}/${queryField}/${queryColumn}/${queryValue}`;
    }
    if ("summarize" in query) {
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

function _urlEncodeWhereQuery(q: WhereQueryTerm): string {
  const { field: queryField, column: queryColumn, value: queryValue } = q;
  return `${_dubEncURIComp(queryField)}:${_dubEncURIComp(
    queryColumn
  )}=${_dubEncURIComp(queryValue)}`;
}

function _urlEncodeSummarizeQuery(q: SummarizeQueryTerm): string {
  const { method, field, column, values } = q;
  const filter = values
    .map((value: string) => _urlEncodeWhereQuery({ field, column, value }))
    .join("&");
  return `method=${method}&${filter}`;
}

export function _urlEncodeComplexQuery(q: ComplexQuery): string {
  if (typeof q === "object") {
    if ("where" in q) {
      return _urlEncodeWhereQuery(q.where);
    }
    if ("summarize" in q) {
      return _urlEncodeSummarizeQuery(q.summarize);
    }
  }
  throw new Error("Unrecognized complex query type");
}

export function _urlEncodeLabelQuery(colKey: string, q: Query): string {
  if (!colKey) throw new Error("Unsupported query by name");
  if (typeof q !== "string") throw new Error("Query must be a simple label.");
  return `${colKey}=${encodeURIComponent(q)}`;
}

/**
 * Generate the column key the server will send us for this query.
 */
export function _hashStringValues(arrayOfString: string[]): string {
  return sha1(arrayOfString.join(""));
}
