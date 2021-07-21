import argparse
import collections
import json
import logging
import math
import string

import anndata
import numpy as np
import pandas as pd
import yaml

from . import gene_symbol
from . import ontology
from . import validate

REPLACE_SUFFIX = "_original"
ONTOLOGY_SUFFIX = "_ontology_term_id"


def is_curie(value):
    """Return True iff the value is an OBO-id CURIE like EFO:000001"""
    return (
        value.count(":")
        and all(len(part) > 0 for part in value.split(":"))
        and all(c in string.digits for c in value.split(":")[1])
    )


def is_ontology_field(field_name):
    """Return True iff the field_name is an ontology field like tissue_ontology_term_id"""
    return field_name.endswith(ONTOLOGY_SUFFIX)


def get_label_field_name(field_name):
    """Get the associated label field from an ontology field, assay_ontology_term_id --> assay"""
    return field_name[: -len(ONTOLOGY_SUFFIX)]


def split_suffix(maybe_curie):
    """Split off the (cell culture) or (organoid) suffix."""

    suffixes = [" (cell culture)", " (organoid)"]
    for suffix in suffixes:
        if maybe_curie.endswith(suffix):
            return maybe_curie[: -len(suffix)], suffix
    return maybe_curie, ""


def get_curie_and_label(maybe_curie):
    """Given a string that might be a curie, return a (curie, label) pair"""

    maybe_curie, suffix = split_suffix(maybe_curie)
    if not is_curie(maybe_curie):
        return ("", maybe_curie + suffix)
    return (maybe_curie + suffix, ontology.get_ontology_label(maybe_curie) + suffix)


def safe_add_field(adata_attr, field_name, field_value):
    """Add a field and value to an AnnData, but don't clobber an exising value."""

    if isinstance(field_value, list) and field_value and isinstance(field_value[0], dict):
        field_value = json.dumps(field_value)
    if field_name in adata_attr:
        adata_attr[field_name + REPLACE_SUFFIX] = adata_attr[field_name]
    adata_attr[field_name] = field_value


def remix_uns(adata, uns_config):
    """Add fields from the config to adata.uns"""
    for field_name, field_value in uns_config.items():

        if is_ontology_field(field_name):
            # If it's an ontology field, look it up
            label_field_name = get_label_field_name(field_name)
            ontology_term, ontology_label = get_curie_and_label(field_value)
            safe_add_field(adata.uns, field_name, ontology_term)
            safe_add_field(adata.uns, label_field_name, ontology_label)
        else:
            safe_add_field(adata.uns, field_name, field_value)


def remix_obs(adata, obs_config):
    """Add fields from the config to adata.obs"""

    for field_name, field_value in obs_config.items():

        if isinstance(field_value, dict):
            # If the value is a dict, that means we are supposed to map from an
            # existing column to the new one
            source_column, column_map = next(iter(field_value.items()))
            nan_value = None
            for key in column_map:
                if isinstance(key, float) and math.isnan(key):
                    nan_value = column_map[key]
            if nan_value is not None:
                column_map["nan"] = nan_value

            for key in column_map:
                if key not in adata.obs[source_column].unique():
                    logging.warning(f'Key {key} not in adata.obs["{source_column}"]')

            for value in adata.obs[source_column].unique():
                if value not in column_map:
                    logging.warning(f'Value {value} in adata.obs["{source_column}"] not in translation dict')

            if is_ontology_field(field_name):
                ontology_term_map, ontology_label_map = {}, {}
                logging.info(f"Looking up labels for {field_name}")
                for original_value, maybe_curie in column_map.items():
                    curie, label = get_curie_and_label(maybe_curie)
                    ontology_term_map[original_value] = curie
                    ontology_label_map[original_value] = label
                    logging.info(f"Mapping {original_value} -> {curie} -> {label}")

                ontology_column = adata.obs[source_column].replace(ontology_term_map, inplace=False)
                label_column = adata.obs[source_column].replace(ontology_label_map, inplace=False)

                safe_add_field(adata.obs, field_name, ontology_column)
                safe_add_field(adata.obs, get_label_field_name(field_name), label_column)
            else:
                label_column = adata.obs[source_column].replace(column_map, inplace=False)
                safe_add_field(adata.obs, field_name, label_column)

        else:
            if is_ontology_field(field_name):
                # If it's an ontology field, look it up
                label_field_name = get_label_field_name(field_name)
                ontology_term, ontology_label = get_curie_and_label(field_value)
                safe_add_field(adata.obs, field_name, ontology_term)
                safe_add_field(adata.obs, label_field_name, ontology_label)
            else:
                safe_add_field(adata.obs, field_name, field_value)


def merge_df(df, domain, index, columns):
    """
    Given a dataframe with duplicate column labels, merge and return a dataframe where
    the duplicates have been merged together, resulting in a dataframe with unique column
    labels.

    "merge" depends on the value of domain. If the domain is "raw", then duplicate columns
    can just be summed. If it's "log1p" or "sqrt", it needs to be exp1m'd or squared, then
    summed, and then logged or sqrt'd again.
    """

    if not isinstance(df, np.ndarray):
        to_merge = df.toarray()
    else:
        to_merge = df
    if domain == "raw":
        merged_df = pd.DataFrame(to_merge, index=index, columns=columns).sum(axis=1, level=0, skipna=False)
    elif domain == "log1p":
        merged_df = pd.DataFrame(np.expm1(to_merge, dtype=np.float128), index=index, columns=columns).sum(
            axis=1, level=0, skipna=False
        )
        merged_df = pd.DataFrame(np.log1p(merged_df.to_numpy()), index=merged_df.index, columns=merged_df.columns)
    elif domain == "sqrt":
        merged_df = pd.DataFrame(np.square(to_merge), index=index, columns=columns).sum(axis=1, level=0, skipna=False)
        merged_df = pd.DataFrame(np.sqrt(merged_df.to_numpy()), index=merged_df.index, columns=merged_df.columns)

    return merged_df


def fixup_gene_symbols(adata, fixup_config):
    """Update the var index to hold a consistent set of HGNC gene symbols."""

    upgraded_var_index = gene_symbol.get_upgraded_var_index(adata.var)

    merged_X = merge_df(adata.X, fixup_config["X"], adata.obs.index, upgraded_var_index)
    fixup_adata = anndata.AnnData(
        X=merged_X,
        obs=adata.obs,
        var=merged_X.columns.to_frame(name="hgnc_gene_symbol"),
        uns=adata.uns,
        obsm=adata.obsm,
    )

    for layer, domain in fixup_config.items():
        if layer == "X":
            continue
        if layer == "raw.X":
            df = adata.raw.X
        else:
            df = adata.layers[layer]

        merged_df = merge_df(df, domain, adata.obs.index, upgraded_var_index)
        assert merged_df.index.equals(merged_X.index)
        assert merged_df.columns.equals(merged_X.columns)

        if domain == "raw":
            fixup_raw = anndata.AnnData(
                X=merged_df,
                obs=adata.obs,
                var=merged_X.columns.to_frame(name="hgnc_gene_symbol"),
            )
            fixup_adata.raw = fixup_raw
        else:
            fixup_adata.layers[layer] = merged_df

    return fixup_adata


def _strip_version(adata):
    """Remove version information from the AnnData object."""

    if "version" in adata.uns_keys():
        del adata.uns["version"]


def apply_schema(source_h5ad, remix_config, output_filename):

    try:
        import scanpy
    except ImportError:
        raise ImportError("scanpy must be installed for cellxgene schema")
    adata = scanpy.read_h5ad(source_h5ad)
    config = yaml.load(open(remix_config), Loader=yaml.FullLoader)
    remix_uns(adata, config["uns"])
    remix_obs(adata, config["obs"])

    if config.get("fixup_gene_symbols"):
        adata = fixup_gene_symbols(adata, config["fixup_gene_symbols"])

    if (
        "version" in adata.uns_keys()
        and isinstance(adata.uns["version"], collections.Mapping)
        and "corpora_schema_version" in adata.uns["version"]
    ):
        schema_version = adata.uns["version"]["corpora_schema_version"]
        try:
            validate.get_schema_definition(schema_version)
        except ValueError:
            logging.warning(
                f"Stripping version information out of AnnData because schema " f"version {schema_version} is unknown."
            )
            _strip_version(adata)

        if not validate.validate_adata(adata, shallow=False):
            logging.warning(
                f"Stripping version information out of AnnData because it does not "
                f"follow schema version {schema_version} ."
            )
            _strip_version(adata)

    adata.write_h5ad(output_filename, compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-h5ad", required=True)
    parser.add_argument("--remix-config", required=True)
    parser.add_argument("--output-filename", required=True)
    args = parser.parse_args()
    apply_schema(args.source_h5ad, args.remix_config, args.output_filename)
