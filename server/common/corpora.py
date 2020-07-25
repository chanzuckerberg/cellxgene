"""
Corpora schema conventions support.  Helper functions for reading.

https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema.md

https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema_h5ad_implementation.md
"""
import json

from server.cli.upgrade import validate_version_str


def corpora_get_versions_from_anndata(adata):
    """
    Given an AnnData object, return:
        * None - if not a Corpora object
        * [ corpora_schema_version, corpora_encoding_version ] - if a Corpora object

    Implements the identification protocol defined in the spect.
    """

    # per Corpora AnnData spec, this is a corpora file if the following is true
    if "version" not in adata.uns or "corpora_schema_version" not in adata.uns["version"]:
        # oops, not corpora
        return None

    corpora_schema_version = adata.uns["version"].get("corpora_schema_version", None)
    corpora_encoding_version = adata.uns["version"].get("corpora_encoding_version", None)

    # TODO: spec says these must be SEMVER values, so check.
    if validate_version_str(corpora_schema_version) and validate_version_str(corpora_encoding_version):
        return [corpora_schema_version, corpora_encoding_version]

    return None


def corpora_get_props_from_anndata(adata):
    """
    Get Corpora dataset properties from an AnnData
    """
    versions = corpora_get_versions_from_anndata(adata)
    if versions is None:
        return None
    [corpora_schema_version, corpora_encoding_version] = versions
    version_is_supported = (
        corpora_schema_version is not None
        and corpora_encoding_version is not None
        and corpora_schema_version.startswith("1.")
        and corpora_encoding_version.startswith("0.1.")
    )

    if not version_is_supported:
        raise ValueError("Unsupported Corpora schema version")

    Required_Simple_Fields = [
        "version",
        "title",
        "layer_descriptions",
        "organism",
        "organism_ontology_term_id",
        "project_name",
        "project_description",
    ]
    # Spec says some values encoded as JSON due to the inability of AnnData to store complex types.
    Required_Json_Fields = ["contributors", "project_links"]
    Optional_Simple_Fields = ["preprint_doi", "publication_doi", "default_embedding", "default_field", "tags"]

    corpora_props = {}
    for key in Required_Simple_Fields:
        if key not in adata.uns:
            raise KeyError(f"missing Corpora schema field {key}")
        corpora_props.update({key: adata.uns[key]})
    for key in Required_Json_Fields:
        if key not in adata.uns:
            raise KeyError(f"missing Corpora schema field {key}")
        corpora_props.update({key: json.loads(adata.uns[key])})
    for key in Optional_Simple_Fields:
        if key in adata.uns:
            corpora_props.update({key: adata.uns[key]})

    return corpora_props
