import React from "react";
import { Colors } from "@blueprintjs/core";

import { AnnotationsHelpers } from "../../util/stateManager";

export function isLabelErroneous(label, metadataField, schema) {
  /*
    return false if this is a LEGAL/acceptable category name or NULL/empty string,
    or return an error type.
    */

  /* allow empty string */
  if (label === "") return false;

  /* check for label syntax errors */
  const error = AnnotationsHelpers.annotationNameIsErroneous(label);
  if (error) return error;

  /* disallow duplicates */
  const { obsByName } = schema.annotations;
  if (obsByName[metadataField].categories.indexOf(label) !== -1)
    return "duplicate";

  /* otherwise, no error */
  return false;
}

/* all other errors - map code to human error message */
const errorMessageMap = {
  "empty-string": "Blank names not allowed",
  duplicate: "Name must be unique",
  "trim-spaces": "Leading and trailing spaces not allowed",
  "illegal-characters":
    "Only alphanumeric and special characters (-_.) allowed",
  "multi-space-run": "Multiple consecutive spaces not allowed",
};

export function labelPrompt(err, prolog, epilog) {
  let errPrompt = null;
  if (err) {
    let errMsg = errorMessageMap[err] ?? "error";
    errMsg = errMsg[0].toLowerCase() + errMsg.slice(1);
    errPrompt = (
      <span
        style={{
          marginTop: 7,
          color: Colors.ORANGE3,
        }}
      >
        {errMsg}
      </span>
    );
  }
  return (
    <span>
      {prolog}
      {err ? " - " : null}
      {errPrompt}
      {epilog}
    </span>
  );
}
