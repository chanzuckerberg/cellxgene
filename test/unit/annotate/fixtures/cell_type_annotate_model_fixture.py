import pickle
import random
from os import PathLike
from tempfile import mkstemp
from typing import Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix


def write_model(model) -> PathLike:
    fn = mkstemp()[1]
    with open(fn, 'wb') as f:
        pickle.dump(model, f)
        return fn


class FakeModel:
    def __init__(self, predicted_labels_, ref_dataset: pd.DataFrame = None):
        self.predicted_labels = predicted_labels_
        self.original_label_key = 'training_label'
        self.adata = ref_dataset

    def predict(self, ad: AnnData) -> Sequence:
        return random.choices(list(self.predicted_labels), k=ad.shape[0])


def write_query_dataset(obs_size: int, gene_identifiers=None) -> PathLike:
    query_dataset_file_path = mkstemp()[1]
    build_dataset(obs_size, gene_identifiers).write_h5ad(query_dataset_file_path)
    return query_dataset_file_path


def build_dataset(obs_size, gene_identifiers=None):
    if gene_identifiers is None:
        gene_identifiers = ['a']
    counts = csr_matrix(np.random.poisson(1, size=(obs_size, len(gene_identifiers))), dtype=np.float32)
    return AnnData(X=counts, var=pd.DataFrame(index=gene_identifiers))
