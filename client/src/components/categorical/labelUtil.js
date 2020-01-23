import React from "react";
import { AnnotationsHelpers } from "../../util/stateManager";

export function isLabelErroneous(label, metadataField, ontology, schema) {
	/*
    return false if this is a LEGAL/acceptable category name or NULL/empty string,
    or return an error type.
    */

	/* allow empty string */
	if (label === "") return false;

	/* check for label syntax errors, but allow terms in ontology */
	const termInOntology = ontology?.termSet.has(label) ?? false;
	const error = AnnotationsHelpers.annotationNameIsErroneous(label);
	if (error && !termInOntology) return error;

	/* disallow duplicates */
	const { obsByName } = schema.annotations;
	if (obsByName[metadataField].categories.indexOf(label) !== -1)
		return "duplicate";

	/* otherwise, no error */
	return false;
}

export function labelErrorMessage(label, metadataField, ontology, schema) {
	const err = isLabelErroneous(label, metadataField, ontology, schema);

	if (err === "duplicate") {
		/* duplicate error is special cased because it has special formatting */
		return (
			<span>
				<span style={{ fontStyle: "italic" }}>{label}</span> already
				exists already exists within{" "}
				<span style={{ fontStyle: "italic" }}>{metadataField}</span>{" "}
			</span>
		);
	}

	if (err) {
		/* all other errors - map code to human error message */
		const errorMessageMap = {
			"empty-string": "Blank names not allowed",
			duplicate: "Name must be unique",
			"trim-spaces": "Leading and trailing spaces not allowed",
			"illegal-characters":
				"Only alphanumeric and special characters (-_.) allowed",
			"multi-space-run": "Multiple consecutive spaces not allowed"
		};
		const errorMessage = errorMessageMap[err] ?? "error";
		return <span>{errorMessage}</span>;
	}

	/* no error, no message generated */
	return null;
}
