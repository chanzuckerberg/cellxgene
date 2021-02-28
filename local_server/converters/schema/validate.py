import json
import re
import os
import sys

import pandas as pd
import yaml


def _is_null(v):
    """Return True if v is null, for one of the multiple ways a "null" value shows up in an h5ad."""
    return pd.isnull(v) or (hasattr(v, "__len__") and len(v) == 0)


def _validate_stringified_list_of_dicts(s):
    """Verify that a string can be parsed into a list.

    We have some types that are lists of dicts. Those cannot be stored directly in an h5ad, so we have to
    json.dumps them. This verifies that we can load them back.
    """

    try:
        list_ = json.loads(s)
        if not isinstance(list_, list):
            return False
        for el in list_:
            if not isinstance(el, dict):
                return False
        return True
    except (json.JSONDecodeError, TypeError):
        pass
    return False


def _validate_human_readable_string(s):
    """Verify that a string is human-readable.

    There are parts of the schema where a "human-readable" string is required. "Human-readable" is kind
    of vague and subjective. I feel like I can read many strings. So here we just check for the main ways
    that fails: someone puts in an ontology term id or and ensembl gene/transcript id.

    Returns False if s is not a string or is one of those bad string types.
    """

    return isinstance(s, str) and (not re.match(r"[A-Z]\w+:\d+", s)) and (not re.match(r"ENS[GT]\d+$", s))


def _validate_curie(c, prefixes):
    """Verify that a string is a valid compact URI, like EFO:000001. If prefixes is not empty, make sure the
    prefix of the curies is in prefixes.
    """

    if not c:
        return True

    match = re.match(r"([A-Z]\w+):\d+$", c)

    if prefixes:
        return match and match.group(1) in prefixes
    else:
        return match


def _validate_suffixed_curie(c, prefixes):
    """Verify that a string is a compact URI with an optional suffix like 'EFO:00001 (cell culture)'"""

    # Pull off the suffix
    suffix = re.findall(r"\ \(.*\)$", c)
    if suffix:
        c = c[: -len(suffix[0])]
    return _validate_curie(c, prefixes)


def _validate_column(column, column_name, df_name, schema_def):
    """Given a schema definition and the column of a dataframe, verify that the column satifies
    the schema.
    """

    errors = []

    if schema_def.get("unique"):
        if column.nunique() != len(column):
            errors.append(f"Column {column_name} in dataframe {df_name} is not unique.")

    if "nullable" in schema_def and not schema_def["nullable"]:
        if any(_is_null(v) for v in column):
            errors.append(f"Column {column_name} in dataframe {df_name} contains empty values.")

    if schema_def.get("type") == "human-readable string":
        non_readables = [v for v in column if not _validate_human_readable_string(v)]
        if non_readables:
            errors.append(
                f"Column {column_name} in dataframe {df_name} contains non-human-readable "
                f"values like {non_readables[0]}"
            )

    if schema_def.get("type") in ("curie", "suffixed curie"):
        validation_func = _validate_curie if schema_def.get("type") == "curie" else _validate_suffixed_curie
        non_valid_curies = [v for v in column if not validation_func(v, schema_def.get("prefixes"))]
        if non_valid_curies:
            errors.append(
                f"Column {column_name} in dataframe {df_name} contains invalid ontology values like "
                f"{non_valid_curies[0]}."
            )
            if "prefixes" in schema_def:
                errors[-1] += f" Values must be curies from one of these ontologies {schema_def['prefixes']}."

    if "enum" in schema_def:
        bad_enums = [v for v in column if v not in schema_def["enum"]]
        if bad_enums:
            errors.append(
                f"Column {column_name} in dataframe {df_name} contains unpermitted values like "
                f"{bad_enums[0]}. Values must be one of {schema_def['enum']}."
            )

    return errors


def _validate_dict(dict_, dict_name, schema_def):
    """Given a schema definition and dict, verify that the dict satifies the schema."""

    errors = []

    for key in schema_def.get("keys", []):
        if key not in dict_:
            errors.append(f"{dict_name} is missing key {key}.")
        elif schema_def["keys"][key]:
            if schema_def["keys"][key]["type"] == "stringified list of dicts":
                if not _validate_stringified_list_of_dicts(dict_[key]):
                    errors.append(
                        f"Key {key} in {dict_name} should be a JSON-encoded list of dicts, but it is {dict_[key]}"
                    )
            elif schema_def["keys"][key]["type"] == "dict":
                errors.extend(_validate_dict(dict_[key], key, schema_def["keys"][key]))
            elif schema_def["keys"][key]["type"] == "curie":
                if not _validate_curie(dict_[key], schema_def["keys"][key]["prefixes"]):
                    errors.append(f"Key {key} in {dict_name} contains invalid ontology value.")
            if "nullable" in schema_def["keys"][key] and not schema_def["keys"][key]["nullable"]:
                if _is_null(dict_[key]):
                    errors.append(f"Key {key} in dict {dict_name} is an empty value.")

    return errors


def _validate_dataframe(df, df_name, schema_def):
    """Given a dataframe and schema definition, verify that the dataframe follows the schema."""

    errors = []

    if "index" in schema_def:
        errors.extend(_validate_column(df.index, "index", df_name, schema_def["index"]))

    for column in schema_def.get("columns", []):
        if column not in df.columns:
            errors.append(f"Dataframe {df_name} is missing column {column}.")
        else:
            errors.extend(_validate_column(df[column], column, df_name, schema_def["columns"][column]))

    return errors


def get_schema_definition(version):
    """Look up and read a schema definition based on a version number like "1.0.0"."""

    path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "schema_definitions", version.replace(".", "_") + ".yaml"
    )

    if not os.path.isfile(path):
        raise ValueError(f"No definition for version {version} found.")

    return yaml.load(open(path), Loader=yaml.FullLoader)


def deep_check(adata, schema_def):
    """Perform a "deep" check of the AnnData object using the schema definition.

    This checks all the columns and unstructured metadata rather than just the version.

    Returns a list of error messages. If that list is empty, the object passed validation.
    """

    errors = []

    for component, component_def in schema_def["components"].items():
        if component_def["type"] == "dataframe":
            errors.extend(_validate_dataframe(getattr(adata, component), component, component_def))
        elif component_def["type"] == "dict":
            errors.extend(_validate_dict(getattr(adata, component), component, component_def))
        else:
            raise ValueError(f"Unexpected component type {component['type']}")

    return errors


def validate_adata(adata, shallow):
    """Validate an AnnData object. If shallow, just check that the required version information is
    present.
    """

    # Does it have the version information written into uns?
    if "version" not in adata.uns_keys() or "corpora_schema_version" not in adata.uns["version"]:
        print("AnnData file is missing corpora version information")
        return False

    # We can stop here if it's a "shallow" check, that is, if we're just
    # checking that version is present.
    if shallow:
        return True

    schema_def = get_schema_definition(adata.uns["version"]["corpora_schema_version"])

    errors = deep_check(adata, schema_def)

    for error in errors:
        print(error)

    return not errors


def validate(h5ad_path, shallow=False):
    """Entry point for validation."""

    try:
        import scanpy
    except ImportError:
        raise ImportError("scanpy must be installed for cellxgene schema")

    try:
        adata = scanpy.read_h5ad(h5ad_path, backed="r")
    except (OSError, TypeError):
        print(f"Unable to open {h5ad_path} with scanpy.")
        sys.exit(1)

    if not validate_adata(adata, shallow):
        sys.exit(1)
