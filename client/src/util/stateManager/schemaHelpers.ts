/*
Helpers for schema management

TODO: all this would be much more natural if done with a framework
like immutable.js
*/
import cloneDeep from "lodash.clonedeep";

import fromEntries from "../fromEntries";
import catLabelSort from "../catLabelSort";

/*
System wide schema assumptions:
  - schema and data wil be consistent (eg, for user-created annotations)
  - schema will be internally self-consistent (eg, index matches columns)
*/

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function indexEntireSchema(schema: any) {
  /* Index schema for ease of use */
  schema.annotations.obsByName = fromEntries(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    schema.annotations?.obs?.columns?.map((v: any) => [v.name, v]) ?? []
  );
  schema.annotations.varByName = fromEntries(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    schema.annotations?.var?.columns?.map((v: any) => [v.name, v]) ?? []
  );
  schema.layout.obsByName = fromEntries(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    schema.layout?.obs?.map((v: any) => [v.name, v]) ?? []
  );
  schema.layout.varByName = fromEntries(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    schema.layout?.var?.map((v: any) => [v.name, v]) ?? []
  );

  return schema;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function _copyObsAnno(schema: any) {
  /* redux copy conventions - WARNING, only for modifying obs annotations */
  return {
    ...schema,
    annotations: {
      ...schema.annotations,
      obs: cloneDeep(schema.annotations.obs),
    },
  };
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function _copyObsLayout(schema: any) {
  return {
    ...schema,
    layout: {
      ...schema.layout,
      obs: cloneDeep(schema.layout.obs),
    },
  };
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function _reindexObsAnno(schema: any) {
  /* reindex obs annotations ONLY */
  schema.annotations.obsByName = fromEntries(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    schema.annotations.obs.columns.map((v: any) => [v.name, v])
  );
  return schema;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function _reindexObsLayout(schema: any) {
  schema.layout.obsByName = fromEntries(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    schema.layout.obs.map((v: any) => [v.name, v])
  );
  return schema;
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function removeObsAnnoColumn(schema: any, name: any) {
  const newSchema = _copyObsAnno(schema);
  newSchema.annotations.obs.columns = schema.annotations.obs.columns.filter(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (v: any) => v.name !== name
  );
  return _reindexObsAnno(newSchema);
}

// @ts-expect-error ts-migrate(6133) FIXME: 'name' is declared but its value is never read.
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function addObsAnnoColumn(schema: any, name: any, defn: any) {
  const newSchema = _copyObsAnno(schema);
  newSchema.annotations.obs.columns.push(defn);
  return _reindexObsAnno(newSchema);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function removeObsAnnoCategory(schema: any, name: any, category: any) {
  /* remove a category from a categorical annotation */
  const categories = schema.annotations.obsByName[name]?.categories;
  if (!categories)
    throw new Error("column does not exist or is not categorical");

  const idx = categories.indexOf(category);
  if (idx === -1) throw new Error("category does not exist");

  const newSchema = _reindexObsAnno(_copyObsAnno(schema));

  /* remove category.  Do not need to resort as this can't change presentation order */
  newSchema.annotations.obsByName[name].categories.splice(idx, 1);
  return newSchema;
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function addObsAnnoCategory(schema: any, name: any, category: any) {
  /* add a category to a categorical annotation */
  const categories = schema.annotations.obsByName[name]?.categories;
  if (!categories)
    throw new Error("column does not exist or is not categorical");

  const idx = categories.indexOf(category);
  if (idx !== -1) throw new Error("category already exists");

  const newSchema = _reindexObsAnno(_copyObsAnno(schema));

  /* add category, retaining presentation sort order */
  const catAnno = newSchema.annotations.obsByName[name];
  catAnno.categories = catLabelSort(catAnno.writable, [
    ...catAnno.categories,
    category,
  ]);
  return newSchema;
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function addObsLayout(schema: any, layout: any) {
  /* add or replace a layout */
  const newSchema = _copyObsLayout(schema);
  newSchema.layout.obs.push(layout);
  return _reindexObsLayout(newSchema);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function removeObsLayout(schema: any, name: any) {
  /* remove a layout */
  const newSchema = _copyObsLayout(schema);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  newSchema.layout.obs = schema.layout.obs.filter((v: any) => v.name !== name);
  return _reindexObsLayout(newSchema);
}
