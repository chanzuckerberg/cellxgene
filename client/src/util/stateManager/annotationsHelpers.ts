/*
Helper functions for user-editable annotations state management.
See also reducers/annotations.js
*/

/*
There are a number of state constraints assumed throughout the
application:
	- all obs annotations are in annoMatrix,
		regardless of whether or not they are user editable.
	- the annoMatrix.schema is always up to date and matches
		the data
	- the schema flag `writable` correctly indicates whether
	  the annotation is editable/mutable.

In addition, the current state management only allows for
categorical annotations to be writable.
*/

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'schema' implicitly has an 'any' type.
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function isCategoricalAnnotation(schema, name) {
  /* 
  we treat any string, categorical or boolean as a categorical.
  Return true/false/undefined (for unkonwn fields)
  */
  const colSchema = schema.annotations.obsByName[name];
  if (colSchema === undefined) return undefined;
  const { type } = colSchema;
  return type === "string" || type === "boolean" || type === "categorical";
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'schema' implicitly has an 'any' type.
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function isContinuousAnnotation(schema, name) {
  /*
  Return true/false/undefined
  */
  const colSchema = schema.annotations.obsByName[name];
  if (colSchema === undefined) return undefined;
  const { type } = colSchema;
  return !(type === "string" || type === "boolean" || type === "categorical");
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'schema' implicitly has an 'any' type.
function _isUserAnnotation(schema, name) {
  return schema.annotations.obsByName[name]?.writable || false;
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'annoMatrix' implicitly has an 'any' typ... Remove this comment to see the full error message
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function isUserAnnotation(annoMatrix, name) {
  return _isUserAnnotation(annoMatrix.schema, name);
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'df' implicitly has an 'any' type.
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
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

const legalCharacters = /^(\w|[ .()-])+$/;
// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'name' implicitly has an 'any' type.
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
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
