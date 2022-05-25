import pickle
from dataclasses import dataclass
from typing import Optional

import anndata
import numpy as np
import pandas as pd
from anndata import AnnData
# from scvi.models import SCANVI
from pandas import DataFrame
from scvi.model import SCANVI
import scanpy as sc

from server.common.utils.data_locator import DataLocator


@dataclass
class CellTypeTissueModel:
    tissue_type_ontology_term_id: str
    model_url: str
    model: SCANVI


# def extract_tissue_types(dataset: AnnData) -> Set[str]:
#     """
#     Returns distinct set of tissue type ontology term IDs.
#     """
#     pass


# @functools.cache
# def lookup_model(model_repository_base_url: str, tissue_type_ontology_term_id: str) -> CellTypeTissueModel:
#     model_url = os.path.join(model_repository_base_url, tissue_type_ontology_term_id, ".pkl")
#     return CellTypeTissueModel(tissue_type_ontology_term_id, model_url=model_url, model=retrieve_model(model_url))


def retrieve_model(model_url: str) -> SCANVI:
    with DataLocator(uri_or_path=model_url).open("rb") as f:
        model = pickle.load(f)
        if not isinstance(model, SCANVI):
            raise "invalid model object at {model_url}: expected type SCANVI, got {model.type}"
        return model


def annotate(query_dataset: AnnData, model: SCANVI, annotation_column_name: str, gene_col_name: Optional[str] = None,
             min_common_gene_count=100) -> None:
    input_dataset = prepare_query_dataset(query_dataset, model, gene_col_name,
                                          min_common_gene_count=min_common_gene_count)

    # TODO: record model metadata (version) in AnnData
    predictions = model.predict(input_dataset)

    query_dataset.obs[annotation_column_name] = pd.Categorical(predictions)
    # TODO: How do obtain confidence scores?
    # query_dataset.obs[f"{annotation_column_name_prefix}_confidence"] = ???


def extract_raw_counts_adata(ad, raw_counts_layer_name) -> AnnData:
    if raw_counts_layer_name:
        X = ad.layers[raw_counts_layer_name]
        var = ad.var
    elif ad.raw:
        X = ad.raw.X
        var = ad.raw.var
    else:
        X = ad.X
        var = ad.var
    return AnnData(X=X, obs=ad.obs, var=var)


# TODO: this should be part of the model's logic; note that build_raw_counts_adata() has been copied from the
#  training pipeline code
def prepare_query_dataset(query_dataset: AnnData, model, gene_col_name: str,
                          raw_counts_input_layer_name: str = None,
                          min_common_gene_count: int = 1500) -> AnnData:
    raw_counts = extract_raw_counts_adata(query_dataset, raw_counts_input_layer_name)
    prepared = make_congruent_genes(raw_counts, model.adata.var, min_common_gene_count, gene_col_name)

    # Add obs annotation column with same name of the one that was used at training time to hold the labeled values.
    # This needs to exist at prediction time, since it existed at training time.
    prepared.obs[model.original_label_key] = ''

    return prepared


def find_common_genes_mask(query_dataset_var: DataFrame, ref_genes_var: DataFrame, gene_col_name):
    """
    Find common genes between query dataset and reference dataset, given the `var` objects from each.
    Uses the gene identifiers from the query dataset that are found in the specified column name, `gene_col_name`.
    If `gene_col_name` is None, then gene identifiers will be taken from `index`. Returns a mask of the
    query_dataset_var rows indicating which rows are in common.
    """
    query_gene_identifiers = query_dataset_var[gene_col_name] if gene_col_name else query_dataset_var.index
    return query_gene_identifiers.isin(ref_genes_var.index)


# adapted from: https://github.com/LungCellAtlas/mapping_data_to_the_HLCA/blob/8bc80dee2c22945dd078514d23916eb186a878c9/scripts/data_import_and_cleaning.py#L29
def make_congruent_genes(query_dataset: AnnData, ref_genes_var: DataFrame, min_common_gene_count: int,
                         gene_col_name: Optional[str] = None) -> AnnData:
    common_genes_mask = find_common_genes_mask(query_dataset.var, ref_genes_var, gene_col_name)

    n_common_genes = np.count_nonzero(common_genes_mask)
    n_ref_genes = ref_genes_var.shape[0]

    if n_common_genes < min_common_gene_count:
        raise ValueError(
                f"There are only {n_common_genes} genes in common between the query and reference datasets, "
                f"but at least {min_common_gene_count} are required")

    print(f"{n_common_genes} genes detected; {n_ref_genes} genes in reference dataset")

    # make query dataset congruent with reference dataset along the var axis, removing and padding genes as needed

    # remove extraneous genes from query dataset
    query_dataset_subset = query_dataset[:, common_genes_mask].copy()

    # add missing genes to query dataset
    if n_common_genes < n_ref_genes:
        # Pad object with 0-valued gene expression vectors, if needed
        print(f'Not all genes were recovered, filling in zeros for {n_ref_genes - n_common_genes} missing genes...')
        genes_missing = set(ref_genes_var.index).difference(set(query_dataset_subset.var_names))
        query_dataset_congruent = pad_dataset_var(query_dataset_subset, genes_missing)
    else:
        query_dataset_congruent = query_dataset_subset

    # re-order query dataset var axis to match reference dataset var axis ordering
    # TODO: Necessary? Does scANVI depend upon ordering of var columns?
    if gene_col_name:
        query_dataset_congruent.var.set_index(gene_col_name, inplace=True)

    return query_dataset_congruent[:, ref_genes_var.index].copy()


# noinspection PyPep8Naming
def pad_dataset_var(query_dataset_subset, genes_to_add):
    X_padded_genes = pd.DataFrame(data=np.zeros((query_dataset_subset.shape[0], len(genes_to_add))),
                                  index=query_dataset_subset.obs_names, columns=genes_to_add)
    adata_padded_genes = sc.AnnData(X_padded_genes)
    query_dataset_congruent = anndata.concat([query_dataset_subset, adata_padded_genes],
                                             axis=1, join='outer',
                                             index_unique=None, merge='unique')
    return query_dataset_congruent
