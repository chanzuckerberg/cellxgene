import pickle
from dataclasses import dataclass

import numpy as np
import pandas as pd
from anndata import AnnData
# from scvi.models import SCANVI
from scvi.model import SCANVI

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


def annotate(query_dataset: AnnData, model: SCANVI, annotation_column_name: str) -> None:
    input_dataset = prepare_query_dataset(query_dataset, model)

    # TODO: necessary to partition into AnnData-per-tissue, then invoke model on each partitioned AnnData object,
    #  then reconstitute into a single AnnData
    #
    # TODO: record model metadata (version) in AnnData
    predictions = model.predict(input_dataset)

    query_dataset.obs[f"{annotation_column_name}"] = pd.Categorical(predictions)
    # TODO: How do obtain confidence scores?
    # query_dataset.obs[f"{annotation_column_name_prefix}_confidence"] = ???


def build_raw_counts_adata(ad, raw_counts_layer_name) -> AnnData:
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
def prepare_query_dataset(query_dataset, model, raw_counts_input_layer_name=None) -> AnnData:
    input_dataset = build_raw_counts_adata(query_dataset, raw_counts_input_layer_name)

    # Add obs annotation column that was used at training time to hold the labeled values. This needs to exist at
    # prediction time, since it existed at training time
    input_dataset.obs[model.original_label_key] = np.zeros(shape=input_dataset.shape[0], dtype=str)

    return input_dataset

