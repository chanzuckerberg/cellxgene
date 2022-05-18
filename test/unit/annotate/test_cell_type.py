import os.path
import pickle
import unittest

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

from server.annotate.cell_type import annotate


class TestAnnotateCellType(unittest.TestCase):
    def test_annotate(self):
        # TODO: figure out how to make this a fixture w/o committing a large model file to the repo.
        #   Currently, just copying a dummy model here
        with open(os.path.join(os.path.dirname(__file__), 'fixtures/model.pkl'), 'rb') as f:
            model = pickle.load(f)

        obs_size = 100
        var_size = model.adata.shape[1]

        counts = csr_matrix(np.random.poisson(1, size=(obs_size, var_size)), dtype=np.float32)
        query_dataset = AnnData(X=counts)

        annotate(query_dataset, model, annotation_column_name_prefix='predicted_cell_type')

        self.assertEquals(obs_size, query_dataset.obs['predicted_cell_type'].shape[0])
        self.assertTrue(isinstance(query_dataset.obs['predicted_cell_type'].dtype,
                                   pd.core.dtypes.dtypes.CategoricalDtype))
        self.assertEquals(0, len(set(query_dataset.obs['predicted_cell_type'].unique()) -
                                 set(model.adata.obs[model.original_label_key].unique())),
                          "only cell types known at training time are added as predicted labels")

    @unittest.skip('unimplemented')
    def test_align_query_dataset_var_axis(self):
        pass
