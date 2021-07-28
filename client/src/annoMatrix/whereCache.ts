/*
Private support functions.

This implements a query resolver cache, mapping a query onto the column labels
resolved by that query. These labels are then used to manage the acutal data cache,
which stores data by the resolved label.

There are three query forms:
  * primitive (string, number) - which is just reference the column label of same value
  * where query (object) - eg, { where: { field: "var", column: "gene", value: "FOXP2" }}
  * summary query (object) - eg, { summarize: { method: "mean", field: "var", column: "gene", values: ["FOXP2", "GNE", "F5"]}}

These queries all resolve to one or more column labels on a field.  This
cache maintains a record of this, allowing direct access to the data caches
without a server round-trip.

The data structure for where queries, the following query against X as an example:
  { where: { field: "var", column: "column_label_in_var", value: "value_in_var_column" } }
results in the following cached entry:
{
  where: {
    X: {
      var: Map(
        column_label_in_var => Map(
          value_in_var_column => [column_label_in_X, ...]
        )
      )
    }
  },
  summarize: {},
}

And for summarize queries, for the following summary on X:
  { summarize: { method: "mean", field: "var", column: "gene", values: ["G1", "G2"]}}
creates a cache entry of:
{
  where: {},
  summarize: {
    X: {
      mean: {
        var: Map(
          "gene" => Map(
            "G1,G2" => [summary_column_label, ...]
          )
        )
      }
    },
  },
}
*/
import { _getColumnDimensionNames } from "./schema";
import { _hashStringValues } from "./query";

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function _whereCacheGet(
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  whereCache: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  schema: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  field: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  query: any
) {
  /* 
	query will either be an where query (object) or a column name (string).

	Return array of column labels or undefined.
	*/

  if (typeof query === "object") {
    if (query.where) {
      const {
        field: queryField,
        column: queryColumn,
        value: queryValue,
      } = query.where;
      const columnMap = whereCache?.where?.[field]?.[queryField];
      return columnMap?.get(queryColumn)?.get(queryValue) ?? [undefined];
    }
    if (query.summarize) {
      const {
        method,
        field: queryField,
        column: queryColumn,
        values: queryValues,
      } = query.summarize;
      const columnMap = whereCache?.summarize?.[field]?.[method]?.[queryField];
      const queryValueHash = _hashStringValues(queryValues);
      return columnMap?.get(queryColumn)?.get(queryValueHash) ?? [undefined];
    }
    return [undefined];
  }

  return _getColumnDimensionNames(schema, field, query) ?? [undefined];
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'field' implicitly has an 'any' type.
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function _whereCacheCreate(field, query, columnLabels) {
  /*
	Create a new whereCache
	*/
  if (typeof query !== "object") return null;

  if (query.where) {
    const {
      field: queryField,
      column: queryColumn,
      value: queryValue,
    } = query.where;
    return {
      where: {
        [field]: {
          [queryField]: new Map([
            [queryColumn, new Map([[queryValue, columnLabels]])],
          ]),
        },
      },
    };
  }
  if (query.summarize) {
    const {
      method,
      field: queryField,
      column: queryColumn,
      values: queryValues,
    } = query.summarize;
    const queryValueHash = _hashStringValues(queryValues);
    return {
      summarize: {
        [field]: {
          [method]: {
            [queryField]: new Map([
              [queryColumn, new Map([[queryValueHash, columnLabels]])],
            ]),
          },
        },
      },
    };
  }
  // oops, not sure what that query is!
  return {};
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'dst' implicitly has an 'any' type.
function __mergeQueries(dst, src) {
  for (const [queryField, columnMap] of Object.entries(src)) {
    dst[queryField] = dst[queryField] || new Map();
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    for (const [queryColumn, valueMap] of columnMap as any) {
      if (!dst[queryField].has(queryColumn))
        dst[queryField].set(queryColumn, new Map());
      for (const [queryValue, columnLabels] of valueMap) {
        dst[queryField].get(queryColumn).set(queryValue, columnLabels);
      }
    }
  }
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'dst' implicitly has an 'any' type.
function __whereCacheMerge(dst, src) {
  /*
	merge src into dst (modifies dst)
	*/
  if (!src || typeof src !== "object") return dst;

  if (src.where) {
    dst.where = dst.where || {};
    for (const [field, query] of Object.entries(src.where)) {
      dst.where[field] = dst.where[field] || {};
      __mergeQueries(dst.where[field], query);
    }
  }
  if (src.summarize) {
    dst.summarize = dst.summarize || {};
    for (const [field, method] of Object.entries(src.summarize)) {
      dst.summarize[field] = dst.summarize[field] || {};
      // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
      for (const [methodName, query] of Object.entries(method)) {
        dst.summarize[field][methodName] =
          dst.summarize[field][methodName] || {};
        __mergeQueries(dst.summarize[field][methodName], query);
      }
    }
  }
  return dst;
}

// @ts-expect-error ts-migrate(7019) FIXME: Rest parameter 'caches' implicitly has an 'any[]' ... Remove this comment to see the full error message
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function _whereCacheMerge(...caches) {
  return caches.reduce(__whereCacheMerge, {});
}
