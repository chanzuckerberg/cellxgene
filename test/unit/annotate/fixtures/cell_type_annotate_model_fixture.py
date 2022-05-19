import pickle
import random
from os import PathLike
from tempfile import mkstemp
from typing import Set, IO, Any, Sequence

import numpy as np
from anndata import AnnData
from numpy import ndarray
from scipy.sparse import csr_matrix


def write_model(model) -> PathLike:
    fn = mkstemp()[1]
    with open(fn, 'wb') as f:
        pickle.dump(model, f)
        return fn


class FakeModel:
    def __init__(self, predicted_labels_):
        self.predicted_labels = predicted_labels_
        self.original_label_key = 'training_label'

    def predict(self, ad: AnnData) -> Sequence:
        return random.choices(list(self.predicted_labels), k=ad.shape[0])


def build_fake_model(predicted_labels: Set[str]) -> Any:

    return FakeModel(predicted_labels)


def write_query_dataset(obs_size: int, var_size: int) -> PathLike:
    query_dataset_file_path = mkstemp()[1]
    build_query_dataset(obs_size, var_size).write_h5ad(query_dataset_file_path)
    return query_dataset_file_path


def build_query_dataset(obs_size, var_size):
    counts = csr_matrix(np.random.poisson(1, size=(obs_size, var_size)), dtype=np.float32)
    return AnnData(X=counts)
