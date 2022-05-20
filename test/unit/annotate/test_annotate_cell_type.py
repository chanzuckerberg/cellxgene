import os
import unittest

import pandas as pd

from server.annotate import cell_type
from server.cli.annotate import _fetch_model
from test.unit.annotate.fixtures.cell_type_annotate_model_fixture import build_query_dataset, build_fake_model


class TestAnnotateCellType(unittest.TestCase):
    def test__fake_model__annotate(self):
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

    @unittest.skipUnless(os.getenv('TEST_REAL_MODEL_URL'), 'remote model url env var missing')
    # For testing with a real model, if available. Models are large, so we don't commit as a test fixture.
    def test__real_model__annotate(self):
        model = _fetch_model(os.getenv('TEST_REAL_MODEL_URL'))

        obs_size = 100
        query_dataset = build_query_dataset(obs_size, model.adata.shape[1])
        query_dataset.var = model.adata.var

        cell_type.annotate(query_dataset, model, annotation_column_name='predicted_cell_type')

        self.assertEquals(obs_size, query_dataset.obs['predicted_cell_type'].shape[0])
        self.assertTrue(isinstance(query_dataset.obs['predicted_cell_type'].dtype,
                                   pd.core.dtypes.dtypes.CategoricalDtype))
        self.assertEquals(0, len(set(query_dataset.obs['predicted_cell_type'].unique()) -
                                 set(model.adata.obs['cell_type'].cat.categories)),
                          "only cell types known at training time are added as predicted labels")

    @unittest.skip('unimplemented')
    def test_align_query_dataset_var_axis(self):
        pass
