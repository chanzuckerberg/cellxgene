import os.path
import pickle
import unittest

import pandas as pd

from server.annotate import cell_type
from test.unit.annotate.fixtures.cell_type_annotate_model_fixture import build_query_dataset, build_fake_model


class TestAnnotateCellType(unittest.TestCase):
    def test_annotate(self):
        predicted_labels = {'x', 'y', 'z'}
        model = build_fake_model(predicted_labels)

        obs_size = 100

        query_dataset = build_query_dataset(obs_size, 10)

        cell_type.annotate(query_dataset, model, annotation_column_name='predicted_cell_type')

        self.assertEquals(obs_size, query_dataset.obs['predicted_cell_type'].shape[0])
        self.assertTrue(isinstance(query_dataset.obs['predicted_cell_type'].dtype,
                                   pd.core.dtypes.dtypes.CategoricalDtype))
        self.assertEquals(0, len(set(query_dataset.obs['predicted_cell_type'].unique()) - predicted_labels),
                          "only cell types known at training time are added as predicted labels")

    @unittest.skip('unimplemented')
    def test_align_query_dataset_var_axis(self):
        pass
