import os
import unittest

import numpy as np
import pandas as pd
import scvi

from server.annotate import cell_type
from test.unit.annotate.fixtures.cell_type_annotate_model_fixture import build_dataset, FakeModel


@unittest.skip(reason="Requires updates to accommodate new model architecture(s)")
class TestAnnotateCellType(unittest.TestCase):
    def test__annotate__adds_annotation_column(self):
        predicted_labels = {'x', 'y', 'z'}
        ref_dataset = build_dataset(100, gene_identifiers=['a', 'b', 'c'])
        model = FakeModel(predicted_labels, ref_dataset)
        query_dataset = build_dataset(10, gene_identifiers=['a', 'b', 'c'])

        obs_predicted = cell_type.predict(query_dataset, model, annotation_column_name='predicted_cell_type')

        self.assertTrue('predicted_cell_type' in obs_predicted.columns)
        self.assertTrue(isinstance(obs_predicted.predicted_cell_type.dtype,
                                   pd.core.dtypes.dtypes.CategoricalDtype),
                        "annotated column is a categorical data type")
        self.assertEquals(0, len(set(obs_predicted.predicted_cell_type.unique()) - predicted_labels),
                          "only cell types known at training time are added as predicted labels")

    @unittest.skipUnless(os.getenv('TEST_REAL_MODEL_URL'), 'remote model url env var missing')
    # For testing with a real model, if available. Models are large, so we don't commit as a test fixture.
    def test__real_model__annotate(self):
        cell_type_classifier_model_path, reference_embedding_model_path = \
            cell_type.fetch_models(os.getenv('TEST_REAL_MODEL_URL'))

        scanvi_model = scvi.model.SCANVI.load(reference_embedding_model_path)

        query_dataset = build_dataset(obs_size=100, gene_identifiers=scanvi_model.adata.var.index)
        query_dataset.var = scanvi_model.adata.var

        obs_predicted = cell_type.predict(query_dataset,
                                          annotation_column_name='predicted_cell_type',
                                          reference_embedding_model_path=reference_embedding_model_path,
                                          cell_type_classifier_model_path=cell_type_classifier_model_path,
                                          )

        self.assertTrue('predicted_cell_type' in obs_predicted.columns)
        self.assertTrue(isinstance(obs_predicted.predicted_cell_type.dtype,
                                   pd.core.dtypes.dtypes.CategoricalDtype),
                        "annotated column is a categorical data type")
        self.assertEqual(100, np.count_nonzero(obs_predicted.predicted_cell_type),
                        "all cells assigned a predicted cell type")
        self.assertEquals(0, len(set(obs_predicted.predicted_cell_type.unique()) -
                                 set(scanvi_model.adata.obs['cell_type'].cat.categories)),
                          "only cell types known at training time are added as predicted labels")

