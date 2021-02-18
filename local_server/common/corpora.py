"""
Corpora schema conventions support.  Helper functions for reading.

https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema.md

https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema_h5ad_implementation.md
"""
import collections
import json

from local_server.cli.upgrade import validate_version_str
from local_server.common.utils.corpora_constants import CorporaConstants


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

    corpora_props = {}
    for key in CorporaConstants.REQUIRED_SIMPLE_METADATA_FIELDS:
        if key not in adata.uns:
            raise KeyError(f"missing Corpora schema field {key}")
        corpora_props[key] = adata.uns[key]

    for key in CorporaConstants.OPTIONAL_JSON_ENCODED_METADATA_FIELD:
        if key not in adata.uns:
            continue
        try:
            corpora_props[key] = json.loads(adata.uns[key])
        except json.JSONDecodeError:
            raise json.JSONDecodeError(f"Corpora schema field {key} is expected to be a valid JSON string")

    for key in CorporaConstants.OPTIONAL_SIMPLE_METADATA_FIELDS:
        if key in adata.uns:
            corpora_props[key] = adata.uns[key]

    return corpora_props
