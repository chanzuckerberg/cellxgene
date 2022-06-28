import os.path
import tempfile
import warnings
from os import path
from typing import List

import scanpy as sc
import xgboost as xgb
from anndata import AnnData
from pandas import DataFrame
from scipy.sparse import spmatrix

from server.common.utils.data_locator import DataLocator

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=UserWarning)

import numpy as np
import pandas as pd


def annotate(query_dataset, annotation_prefix, model_cache_dir, model_url,
             counts_layer, gene_column_name, use_gpu: bool = False, train_param_overrides: bool = False):
    cell_type_classifier_model_locator = DataLocator(path.join(model_url, 'classifier.xgb'),
                                                     local_cache_dir=model_cache_dir)
    reference_embedding_model_locator = DataLocator(path.join(model_url, 'model.pt'),
                                                    local_cache_dir=model_cache_dir)

    # Retrieve (possibly) remote model file and cache locally. Then place the model file in a temp dir with the
    # specific 'model.pt' file name in order to appease scArches' model loading logic, which expects the directory of
    # the model, not the full path to the model.pt file. While scArches appears to have it own support for caching
    # remote files, it doesn't seem to work.
    with reference_embedding_model_locator.open() as f_in, tempfile.TemporaryDirectory() as tmp_ref_emb_model_dir:
        with open(os.path.join(tmp_ref_emb_model_dir, 'model.pt'), "wb") as f_out:
            f_out.write(f_in.read())

        query_dataset_prepped = prep_query_data(query_dataset,
                                                tmp_ref_emb_model_dir,
                                                counts_layer=counts_layer,
                                                gene_column_name=gene_column_name,
                                                dataset_name=str(query_dataset.filename))

        query_dataset_reference_embedding = map_to_ref(query_dataset_prepped,
                                                       tmp_ref_emb_model_dir,
                                                       use_gpu,
                                                       train_kwargs_overrides=(early_stopping_kwargs_scarches |
                                                                               train_param_overrides))
        query_dataset.obsm[annotation_prefix] = query_dataset_reference_embedding

    obs_predictions = predict(query_dataset_prepped,
                              query_dataset_reference_embedding,
                              cell_type_classifier_model_locator,
                              annotation_column_name=f"{annotation_prefix}_predicted",
                              use_gpu=use_gpu)
    query_dataset.obs = pd.concat([query_dataset.obs, obs_predictions], axis=1)

    generate_query_dataset_umap(query_dataset, annotation_prefix)

    record_prediction_run_metadata(query_dataset, annotation_prefix, model_url)


def generate_query_dataset_umap(query_dataset, annotation_prefix):
    neighbors_key = f"{annotation_prefix}_neighbors"
    sc.pp.neighbors(query_dataset, use_rep=annotation_prefix, key_added=neighbors_key)
    query_dataset_with_umap = sc.tl.umap(query_dataset, neighbors_key=neighbors_key, copy=True)
    query_dataset.obsm[f"X_{annotation_prefix}_umap"] = query_dataset_with_umap.obsm['X_umap']


def record_prediction_run_metadata(query_dataset: AnnData, annotation_prefix: str, model_url: str):
    query_dataset.uns[f"{annotation_prefix}_predictions"] = {'model_url': model_url}


def predict(query_dataset: AnnData,
            query_dataset_embedding,
            cell_type_classifier_model_locator: DataLocator,
            annotation_column_name: str,
            use_gpu: bool = False
            ) -> DataFrame:
    """
    Returns `obs` DataFrame with prediction columns appended
    """
    tree_method = 'gpu_hist' if use_gpu else 'hist'
    cell_type_classifier_model = xgb.XGBClassifier(tree_method=tree_method, objective='multi:softprob')
    with cell_type_classifier_model_locator.local_handle() as fname:
        cell_type_classifier_model.load_model(fname=fname)

    project_labels(query_dataset, cell_type_classifier_model, query_dataset_embedding,
                   annotation_column_name=annotation_column_name)

    return query_dataset.obs


# Parameters
early_stopping_kwargs_scarches = {
    'enable_progress_bar': True,  # Fix for when using torch nightly build
    'early_stopping': True,
    'early_stopping_monitor': 'elbo_train',
    'early_stopping_patience': 10,
    'early_stopping_min_delta': 0.001,
}

plan_kwargs = {
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
    "weight_decay": 0,
}


def extract_raw_counts_adata(ad, raw_counts_layer_name) -> spmatrix:
    if raw_counts_layer_name:
        return ad.layers[raw_counts_layer_name]
    elif ad.raw:
        return ad.raw.X
    else:
        return ad.X


def extract_gene_names(ad, gene_column_name) -> List[str]:
    var_names = ad.var[gene_column_name] if gene_column_name else ad.var_names
    var_names_upcased = [n.upper() for n in var_names]  # TODO: probably wrong to upcase, are HGNC names case-sensitive?
    return var_names_upcased


def build_query_dataset(adata, counts_layer, gene_column_name):
    raw_counts = extract_raw_counts_adata(adata, counts_layer)
    gene_names = extract_gene_names(adata, gene_column_name)
    return AnnData(X=raw_counts,
                   var=pd.DataFrame(index=gene_names),
                   obs=pd.DataFrame(index=adata.obs_names))


def prep_query_data(adata,
                    model_dir,
                    counts_layer='counts',
                    gene_column_name=None,
                    dataset_name='query_data'):  # prep adata object for ref mapping; output -> prepped query object
    """
    A function to prepare the `AnnData` object for reference mapping to the HCLA

    Note: other reference models may have other requirements, such as:
        * count data stored in a different location (not a layer?)
        * size_factors stored or to be computed
        * feature names no HGNC gene names

    Input:
        * `adata`: A query AnnData object with count data stored in adata.layers
        * `model_path`: The HLCA reference model directory
        * `counts_layer`: The name of the layer containing count data. If None, uses raw.X, if available, otherwise X
        * `dataset_name`: (Optional) The name of the dataset

    Output:
        A prepared `AnnData` object that can be loaded into `map_to_ref` for reference mapping
    """
    import scanpy as sc
    import scvi

    # Note: ct_key, batch_key, and unlabeled_category are reference-model-specific
    #       => this would need to be stored somewhere (or standardized)
    ct_key = 'ann_finest_level'
    batch_key = "dataset"
    unlabeled_category = "unlabeled"

    # Prep query anndata
    # Checks:
    # - cells have >= 200 counts -> done
    # - layers['counts'] exists -> done
    # - test type of input

    # To implement:
    # - var names are gene names -> how to test? Just check some examples?
    # - layer contains non-negative integers
    if not isinstance(adata, sc.AnnData):
        raise TypeError(f'`adata` should contain an AnnData object!')

    if not isinstance(dataset_name, str):
        raise TypeError('`dataset_name` should be a string!')

    adata_query = build_query_dataset(adata, counts_layer, gene_column_name)

    if np.any(adata_query.X.sum(1) < 200):  # TODO: parameterize
        raise ValueError('Cells with fewer than 200 counts in the data.\n'
                         'Please filter your data with `sc.pp.filter_cells()` prior to reference mapping.')

    adata_query.obs[batch_key] = dataset_name
    adata_query.obs[ct_key] = unlabeled_category

    # Prepare AnnData object to match to scANVI query data
    scvi.model.SCANVI.prepare_query_anndata(adata_query, model_dir)

    # Add back the required layer
    adata_query.layers['counts'] = adata_query.X

    return adata_query


def map_to_ref(query_dataset,
               model_path,
               use_gpu=False,
               train_kwargs_overrides={},
               plan_kwargs_overrides={},
               ) -> np.ndarray:
    """
    A function to map the query data to the reference atlas

    Input:
        * query_dataset: An AnnData object that has been prepped by the `prep_query_data` function
        * model_path: The HLCA reference model directory
        * max_epochs: The maximum number of epochs for reference training (Default: 500)
        * use_gpu: Boolean flag whether to use gpu for training (Default: True)
        * train_kwargs: Keyword arguments passed to `scvi.model.SCANVI.train()` for early stopping
        * plan_kwargs: Dictionary of training plan keyword arguments passed to `scvi.model.SCANVI.train()`

    Output:
        A numpy.ndarray containing the embedding coordinates of the query data in the HLCA latent space.
        If `output_model` is set to True, then the trained reference model is also output

    """
    import scvi

    # Load query data into the model
    vae_q = scvi.model.SCANVI.load_query_data(
            query_dataset,
            model_path,
            freeze_dropout=True,
    )

    # Train scArches model for query mapping
    vae_q.train(
            max_epochs=500,
            plan_kwargs=(plan_kwargs | plan_kwargs_overrides),
            **(early_stopping_kwargs_scarches | train_kwargs_overrides),
            check_val_every_n_epoch=1,
            use_gpu=use_gpu # 'mps:0'  # TODO: use_gpu,
    )

    emb = vae_q.get_latent_representation()

    return emb


def project_labels(query_dataset: AnnData,
                   cell_type_classifier_model: xgb.XGBClassifier,
                   embedding_coords: np.ndarray,
                   annotation_column_name='pred_labels'):
    """
    A function that projects predicted labels onto the query dataset, along with uncertainty scores.
    Performs in-place update of the adata object, adding columns to the `obs` DataFrame.

    Input:
        * query_dataset: The query `AnnData` object which was mapped to the HLCA reference
        * model_file: Path to the classification model file
        * embedding_coords: A numpy.ndarray of the latent space embedding coordinates of the
                            query data mapped to the HLCA reference. This is output by the
                            `map_to_ref()` function
        * prediction_key: Column name in `adata.obs` where to store the predicted labels

    Output:
        Nothing is output, the passed anndata is modified inplace

    """
    # TODO: uncertainty threshold 0.5

    # Set up query data
    q_emb_df = pd.DataFrame(data=embedding_coords, index=query_dataset.obs_names)

    # Predict labels
    query_dataset.obs[annotation_column_name] = cell_type_classifier_model.predict(q_emb_df)

    # Predict probabilities
    probs = cell_type_classifier_model.predict_proba(q_emb_df)

    # Format probabilities
    df_probs = pd.DataFrame(probs, columns=cell_type_classifier_model.classes_, index=query_dataset.obs_names)

    # TODO: Check this! Does that selected max correspond to the label that was predicted?
    query_dataset.obs[annotation_column_name + "_uncertainty"] = 1 - df_probs.max(1)

    # TODO: Update caller to copy this to output AnnData object, if we desired this
    # adata.uns[annotation_column_name + '_uncertainty'] = df_probs
