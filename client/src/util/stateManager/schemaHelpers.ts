/*
Helpers for schema management

TODO: all this would be much more natural if done with a framework
like immutable.js
*/
import cloneDeep from "lodash.clonedeep";

import fromEntries from "../fromEntries";
import catLabelSort from "../catLabelSort";
import {
  RawSchema,
  Schema,
  LayoutColumn,
  AnnotationColumn,
} from "../../common/types/entities";

/*
System wide schema assumptions:
  - schema and data wil be consistent (eg, for user-created annotations)
  - schema will be internally self-consistent (eg, index matches columns)
*/

export function indexEntireSchema(schema: RawSchema): Schema {
  /* Index schema for ease of use */
  (schema as Schema).annotations.obsByName = fromEntries(
    schema.annotations?.obs?.columns?.map((v) => [v.name, v]) || []
  );
  (schema as Schema).annotations.varByName = fromEntries(
    schema.annotations?.var?.columns?.map((v) => [v.name, v]) || []
  );
  (schema as Schema).layout.obsByName = fromEntries(
    schema.layout?.obs?.map((v) => [v.name, v]) || []
  );
  (schema as Schema).layout.varByName = fromEntries(
    schema.layout?.var?.map((v) => [v.name, v]) || []
  );

  return schema as Schema;
}

function _copyObsAnno(schema: Schema): Schema {
  /* redux copy conventions - WARNING, only for modifying obs annotations */
  return {
    ...schema,
    annotations: {
      ...schema.annotations,
      obs: cloneDeep(schema.annotations.obs),
    },
  };
}

function _copyObsLayout(schema: Schema): Schema {
  return {
    ...schema,
    layout: {
      ...schema.layout,
      obs: cloneDeep(schema.layout.obs),
    },
  };
}

function _reindexObsAnno(schema: Schema): Schema {
  /* reindex obs annotations ONLY */
  schema.annotations.obsByName = fromEntries(
    schema.annotations.obs.columns.map((v) => [v.name, v])
  );
  return schema;
}

function _reindexObsLayout(schema: Schema) {
  schema.layout.obsByName = fromEntries(
    schema.layout.obs.map((v) => [v.name, v])
  );
  return schema;
}

export function removeObsAnnoColumn(schema: Schema, name: string): Schema {
  const newSchema = _copyObsAnno(schema);
  newSchema.annotations.obs.columns = schema.annotations.obs.columns.filter(
    (v) => v.name !== name
  );
  return _reindexObsAnno(newSchema);
}

export function addObsAnnoColumn(
  schema: Schema,
  _: string,
  defn: AnnotationColumn
): Schema {
  const newSchema = _copyObsAnno(schema);

  newSchema.annotations.obs.columns.push(defn);

  return _reindexObsAnno(newSchema);
}

export function removeObsAnnoCategory(
  schema: Schema,
  name: string,
  category: string
): Schema {
  /* remove a category from a categorical annotation */
  const categories = schema.annotations.obsByName[name]?.categories;

  if (!categories) {
    throw new Error("column does not exist or is not categorical");
  }

  const idx = categories.indexOf(category);

  if (idx === -1) throw new Error("category does not exist");

  const newSchema = _reindexObsAnno(_copyObsAnno(schema));

  /* remove category.  Do not need to resort as this can't change presentation order */
  newSchema.annotations.obsByName[name].categories?.splice(idx, 1);

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
    ...(catAnno.categories || []),
    category,
  ]);

  return newSchema;
}

export function addObsLayout(schema: Schema, layout: LayoutColumn): Schema {
  /* add or replace a layout */
  const newSchema = _copyObsLayout(schema);
  newSchema.layout.obs.push(layout);
  return _reindexObsLayout(newSchema);
}

export function removeObsLayout(schema: Schema, name: string): Schema {
  /* remove a layout */
  const newSchema = _copyObsLayout(schema);

  newSchema.layout.obs = schema.layout.obs.filter((v) => v.name !== name);

  return _reindexObsLayout(newSchema);
}
