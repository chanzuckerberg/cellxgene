/*
Private support functions.

Support for a "where" query, eg,

	{ where: { field: "var", column: "gene", value: "FOXP2" }}

These evaluate to a given column label.

The "where cache" is a map that saves evaluated queries and points 
to the column label they resolve to.

Data structure, using X as the example field being queried, and var as
the index.

{ 
	X: {
		var: Map(
			column_label_in_var => Map(value_in_var_column => [column_label_in_X, ...])
		)
	}
}
*/
import { _getColumnDimensionNames } from "./schema";

export function _whereCacheGet(whereCache, schema, field, query) {
  /* 
	query will either be an where query (object) or a column name (string).

	Return array of column labels or undefined.
	*/

  if (typeof query === "object") {
    const { field: queryField, column: queryColumn, value: queryValue } = query;

    const columnMap = whereCache?.[field]?.[queryField];
    if (columnMap === undefined) return [undefined];

    const valueMap = columnMap.get(queryColumn);
    if (valueMap === undefined) return [undefined];

    const columnLabels = valueMap.get(queryValue);
    return columnLabels === undefined ? [undefined] : columnLabels;
  }

  const colDims = _getColumnDimensionNames(schema, field, query);
  return colDims === undefined ? [undefined] : colDims;
}

export function _whereCacheCreate(field, query, columnLabels) {
  /*
	Create a new whereCache
	*/
  if (typeof query !== "object") return null;

  const { field: queryField, column: queryColumn, value: queryValue } = query;
  const whereCache = {
    [field]: {
      [queryField]: new Map([
        [queryColumn, new Map([[queryValue, columnLabels]])],
      ]),
    },
  };
  return whereCache;
}

function __whereCacheMerge(dst, src) {
  /*
	merge src into dst (modifies dst)
	*/
  if (!dst) dst = {};
  if (!src || typeof src !== "object") return dst;
  Object.entries(src).forEach(([field, query]) => {
    if (!Object.prototype.hasOwnProperty.call(dst, field)) dst[field] = {};
    Object.entries(query).forEach(([queryField, columnMap]) => {
      if (!Object.prototype.hasOwnProperty.call(dst[field], queryField))
        dst[field][queryField] = new Map();
      columnMap.forEach((valueMap, queryColumn) => {
        if (!dst[field][queryField].has(queryColumn))
          dst[field][queryField].set(queryColumn, new Map());
        valueMap.forEach((columnLabels, queryValue) => {
          dst[field][queryField].get(queryColumn).set(queryValue, columnLabels);
        });
      });
    });
  });
  return dst;
}

export function _whereCacheMerge(...caches) {
  return caches.reduce((dst, src) => __whereCacheMerge(dst, src), {});
}
