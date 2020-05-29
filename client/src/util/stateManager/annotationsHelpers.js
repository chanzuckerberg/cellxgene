/*
Helper functions for user-editable annotations state management.
See also reducers/annotations.js
*/
import { unassignedCategoryLabel } from "../../globals";
import * as SchemaHelpers from "./schemaHelpers";
import { obsAnnoDimensionName } from "../nameCreators";

/*
There are a number of state constraints assumed throughout the
application:
	- all obs annotations are in {world|universe}.obsAnnotations,
		regardless of whether or not they are user editable.
	- the {world|universe}.schema is always up to date and matches
		the data
	- the schema flag `writable` correctly indicates whether
	  the annotation is editable/mutable.

In addition, the current state management only allows for
categorical annotations to be writable.
*/

export function isCategoricalAnnotation(schema, name) {
  /* we treat any string, categorical or boolean as a categorical */
  const { type } = schema.annotations.obsByName[name];
  return type === "string" || type === "boolean" || type === "categorical";
}

export function isContinuousAnnotation(schema, name) {
  return !isCategoricalAnnotation(schema, name);
}

function _isUserAnnotation(schema, name) {
  return schema.annotations.obsByName[name]?.writable;
}

export function isUserAnnotation(worldOrUniverse, name) {
  return _isUserAnnotation(worldOrUniverse.schema, name);
}

export function removeObsAnnoSchema(schema, name) {
  /*
	remove named annotation from obs annotation schema
	*/

  /* only remove if it exists and is a user annotation */
  if (!_isUserAnnotation(schema, name))
    throw new Error("removing non-user-defined schema");
  return SchemaHelpers.removeObsAnnoColumn(schema, name);
}

export function addObsAnnoSchema(schema, name, colSchema) {
  /*
	add a categorical type to the obs annotation schema 
	*/

  /* collision detection */
  if (schema.annotations.obs.columns.some((v) => v.name === name))
    throw Error("annotations may not contain duplicate category names");
  if (name !== colSchema.name) throw Error("column schema does not match");
  return SchemaHelpers.addObsAnnoColumn(schema, name, colSchema);
}

export function dupObsAnnoSchema(schema, sourceName, dupName, defaultSchema) {
  /*
	duplicate the obs annotation `sourceName` schema, but with the name `dupName`
	*/
  const colSchema = {
    ...schema.annotations.obsByName[sourceName],
    ...defaultSchema,
    name: dupName,
  };
  /* existance check */
  if (!colSchema) throw Error("source annotation does not exist");
  /* collision detection */
  if (schema.annotations.obs.columns.some((v) => v.name === dupName))
    throw Error("annotations may not contain duplicate category names");
  return SchemaHelpers.addObsAnnoColumn(schema, dupName, colSchema);
}

export function removeObsAnnoCategory(schema, name, category) {
  /* don't allow deletion of unassigned category on writable annotations */

  if (!_isUserAnnotation(schema, name))
    throw new Error("unable to modify read-only schema");
  if (category === unassignedCategoryLabel)
    throw new Error("may not remove unassigned category label");

  return SchemaHelpers.removeObsAnnoCategory(schema, name, category);
}

export function addObsAnnoCategory(schema, name, category) {
  if (!_isUserAnnotation(schema, name))
    throw new Error("unable to modify read-only schema");

  return SchemaHelpers.addObsAnnoCategory(schema, name, category);
}

export function setLabelByValue(df, colName, fromLabel, toLabel) {
  /* 
	in the dataframe column `colName`, set any value of `fromLabel` to `toLabel`
	*/
  const keys = df.colIndex.labels();
  const ndf = df.mapColumns((col, colIdx) => {
    if (colName !== keys[colIdx]) return col;

    /* clone data and return it. */
    const newCol = col.slice();
    for (let i = 0, l = newCol.length; i < l; i += 1) {
      if (newCol[i] === fromLabel) newCol[i] = toLabel;
    }
    return newCol;
  });
  return ndf;
}

export function setLabelByMask(df, colName, mask, label) {
  /* 
	in the dataframe column `colName`, set the masked rows to 'label'
	*/
  const keys = df.colIndex.labels();
  const ndf = df.mapColumns((col, colIdx) => {
    if (colName !== keys[colIdx]) return col;

    /* clone data and return it. */
    const newCol = col.slice();
    for (let i = 0, l = newCol.length; i < l; i += 1) {
      if (mask[i]) newCol[i] = label;
    }
    return newCol;
  });
  return ndf;
}

export function allHaveLabelByMask(df, colName, label, mask) {
  // return true if all rows as indicated by mask have the colname set to label.
  // False if not.
  const col = df.col(colName);
  if (!col) return false;
  if (df.length !== mask.length)
    throw new RangeError("mismatch on mask length");

  for (let i = 0; i < df.length; i += 1) {
    if (mask[i]) {
      if (col.iget(i) !== label) return false;
    }
  }
  return true;
}

export function worldToUniverseMask(worldMask, worldObsAnnotations, nObs) {
  /*
	given world seleciton mask, return a selection mask for entire universe 
	that has same selection state.
	*/
  const mask = new Uint8Array(nObs);
  const { rowIndex } = worldObsAnnotations;

  for (let i = 0, l = worldMask.length; i < l; i += 1) {
    if (worldMask[i]) {
      const label = rowIndex.getLabel(i);
      mask[label] = 1;
    }
  }
  return mask;
}

export function createWritableAnnotationDimensions(world, crossfilter) {
  const { obsAnnotations, schema } = world;
  const writableAnnotations = schema.annotations.obs.columns
    .filter((s) => s.writable)
    .map((s) => s.name);

  crossfilter = writableAnnotations.reduce((xflt, anno) => {
    const dimName = obsAnnoDimensionName(anno);
    if (xflt.hasDimension(dimName)) xflt = xflt.delDimension(dimName);
    return xflt.addDimension(
      dimName,
      "enum",
      obsAnnotations.col(anno).asArray()
    );
  }, crossfilter);
  return crossfilter;
}

const legalCharacters = /^(\w|[ .()-])+$/;
export function annotationNameIsErroneous(name) {
  /*
	Validate the name - return:
	* false - a valid name
	* string - a named error, indicating why it was invalid.

	Tests:
	0. must be string, non-null
	1. no leading or trailing spaces
	2. only accept alpha, numeric, underscore, period, parens, hyphen and space
	3. no runs of multiple spaces
	*/
  if (name === "") {
    return "empty-string";
  }
  if (name[0] === " " || name[name.length - 1] === " ") {
    return "trim-spaces";
  }
  if (!legalCharacters.test(name)) {
    return "illegal-characters";
  }
  for (let i = 1, l = name.length; i < l; i += 1) {
    if (name[i] === " " && name[i - 1] === " ") {
      return "multi-space-run";
    }
  }

  /* all is well!  Indicte not erroneous with a false */
  return false;
}
