"""
Corpora schema conventions support.  Helper functions for reading.

https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema.md

https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema_h5ad_implementation.md
"""
import collections
import json

from server.cli.upgrade import validate_version_str


def corpora_get_versions_from_anndata(adata):
    """
    Given an AnnData object, return:
        * None - if not a Corpora object
        * [ corpora_schema_version, corpora_encoding_version ] - if a Corpora object

    Implements the identification protocol defined in the specification.
    """

    # per Corpora AnnData spec, this is a corpora file if the following is true
    if "version" not in adata.uns_keys():
        return None
    version = adata.uns["version"]
    if not isinstance(version, collections.abc.Mapping) or "corpora_schema_version" not in version:
        return None

    corpora_schema_version = version.get("corpora_schema_version")
    corpora_encoding_version = version.get("corpora_encoding_version")

    # TODO: spec says these must be SEMVER values, so check.
    if validate_version_str(corpora_schema_version) and validate_version_str(corpora_encoding_version):
        return [corpora_schema_version, corpora_encoding_version]


def corpora_is_version_supported(corpora_schema_version, corpora_encoding_version):
    return (
        corpora_schema_version
        and corpora_encoding_version
        and corpora_schema_version.startswith("1.")
        and corpora_encoding_version.startswith("0.1.")
    )


def corpora_get_props_from_anndata(adata):
    """
    Get Corpora dataset properties from an AnnData
    """
    versions = corpora_get_versions_from_anndata(adata)
    if versions is None:
        return None
    [corpora_schema_version, corpora_encoding_version] = versions
    version_is_supported = corpora_is_version_supported(corpora_schema_version, corpora_encoding_version)
    if not version_is_supported:
        raise ValueError("Unsupported Corpora schema version")

    required_simple_fields = [
        "version",
        "title",
        "layer_descriptions",
        "organism",
        "organism_ontology_term_id",
        "project_name",
        "project_description",
    ]
    # Spec says some values encoded as JSON due to the inability of AnnData to store complex types.
    required_json_fields = ["contributors", "project_links"]
    optional_simple_fields = ["preprint_doi", "publication_doi", "default_embedding", "default_field", "tags"]

    corpora_props = {}
    for key in required_simple_fields:
        if key not in adata.uns:
            raise KeyError(f"missing Corpora schema field {key}")
        corpora_props[key] = adata.uns[key]

    for key in required_json_fields:
        if key not in adata.uns:
            raise KeyError(f"missing Corpora schema field {key}")
        try:
            corpora_props[key] = json.loads(adata.uns[key])
        except json.JSONDecodeError:
            raise json.JSONDecodeError(f"Corpora schema field {key} is expected to be a valid JSON string")

    for key in optional_simple_fields:
        if key in adata.uns:
            corpora_props[key] = adata.uns[key]

    return corpora_props
