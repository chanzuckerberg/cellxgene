import os
import unittest

import numpy as np
import pandas as pd

from server.annotate import cell_type
from server.cli.annotate import _fetch_model
from test.unit.annotate.fixtures.cell_type_annotate_model_fixture import build_dataset, FakeModel


class TestAnnotateCellType(unittest.TestCase):
    def test__annotate__adds_annotation_column(self):
        predicted_labels = {'x', 'y', 'z'}
        ref_dataset = build_dataset(100, gene_identifiers=['a', 'b', 'c'])
        model = FakeModel(predicted_labels, ref_dataset)
        query_dataset = build_dataset(10, gene_identifiers=['a', 'b', 'c'])

        cell_type.annotate(query_dataset, model, annotation_column_name='predicted_cell_type', min_common_gene_count=3)

        self.assertTrue('predicted_cell_type' in query_dataset.obs.columns)
        self.assertTrue(isinstance(query_dataset.obs['predicted_cell_type'].dtype,
                                   pd.core.dtypes.dtypes.CategoricalDtype),
                        "annotated column is a categorical data type")
        self.assertEquals(0, len(set(query_dataset.obs['predicted_cell_type'].unique()) - predicted_labels),
                          "only cell types known at training time are added as predicted labels")

    @unittest.skipUnless(os.getenv('TEST_REAL_MODEL_URL'), 'remote model url env var missing')
    # For testing with a real model, if available. Models are large, so we don't commit as a test fixture.
    def test__real_model__annotate(self):
        model = _fetch_model(os.getenv('TEST_REAL_MODEL_URL'))
        query_dataset = build_dataset(obs_size=100, gene_identifiers=model.adata.var.index)
        query_dataset.var = model.adata.var

        cell_type.annotate(query_dataset, model, annotation_column_name='predicted_cell_type')

        self.assertTrue('predicted_cell_type' in query_dataset.obs.columns)
        self.assertTrue(isinstance(query_dataset.obs['predicted_cell_type'].dtype,
                                   pd.core.dtypes.dtypes.CategoricalDtype),
                        "annotated column is a categorical data type")
        self.assertEqual(100, np.count_nonzero(query_dataset.obs['predicted_cell_type']),
                        "all cells assigned a predicted cell type")
        self.assertEquals(0, len(set(query_dataset.obs['predicted_cell_type'].unique()) -
                                 set(model.adata.obs['cell_type'].cat.categories)),
                          "only cell types known at training time are added as predicted labels")

    def test__make_congruent_genes__adds_missing_genes(self):
        query_dataset = build_dataset(10, ['a', 'c', 'e'])
        query_dataset.var['gene_name'] = ['gene_a', 'gene_c', 'gene_e']
        ref_genes_var = pd.DataFrame(index=['a', 'b', 'c', 'd', 'e'])

        query_dataset_congruent_genes = cell_type.make_congruent_genes(query_dataset, ref_genes_var,
                                                                       min_common_gene_count=3)

        self.assertListEqual(['a', 'b', 'c', 'd', 'e'], query_dataset_congruent_genes.var.index.to_list())

    def test__make_congruent_genes__removes_extraneous_genes(self):
        query_dataset = build_dataset(10, ['a', 'b', 'c', 'd', 'e'])
        query_dataset.var['gene_name'] = ['gene_a', 'gene_b', 'gene_c', 'gene_d', 'gene_e']
        ref_genes_var = pd.DataFrame(index=['a', 'c', 'e'])

        query_dataset_congruent_genes = cell_type.make_congruent_genes(query_dataset, ref_genes_var,
                                                                       min_common_gene_count=3)

        self.assertListEqual(['a', 'c', 'e'], query_dataset_congruent_genes.var.index.to_list())

    def test__make_congruent_genes__misordered_genes_are_ordered(self):
        query_dataset = build_dataset(10, ['a', 'b', 'c', 'd', 'e'])
        query_dataset.var['gene_name'] = ['gene_a', 'gene_d', 'gene_c', 'gene_b', 'gene_e']
        ref_genes_var = pd.DataFrame(index=['a', 'b', 'c', 'd', 'e'])

        query_dataset_congruent_genes = cell_type.make_congruent_genes(query_dataset, ref_genes_var,
                                                                       min_common_gene_count=5)

        self.assertListEqual(['a', 'b', 'c', 'd', 'e'], query_dataset_congruent_genes.var.index.to_list())

    def test__make_congruent_genes__error_on_insufficient_common_genes(self):
        query_dataset = build_dataset(10, ['a', 'b', 'x', 'y', 'z'])
        ref_genes_var = pd.DataFrame(index=['a', 'b', 'c', 'd', 'e'])

        with self.assertRaises(ValueError):
            cell_type.make_congruent_genes(query_dataset, ref_genes_var, min_common_gene_count=3)

    def test__make_congruent_genes__uses_specified_gene_col_name(self):
        query_dataset = build_dataset(10, [0, 1, 2, 3, 4])
        query_dataset.var['gene_name'] = ['a', 'b', 'c', 'd', 'e']
        ref_genes_var = pd.DataFrame(index=['a', 'c', 'e'])

        query_dataset_congruent_genes = cell_type.make_congruent_genes(query_dataset, ref_genes_var,
                                                                       min_common_gene_count=3,
                                                                       gene_col_name='gene_name')

        self.assertListEqual(['a', 'c', 'e'], query_dataset_congruent_genes.var.index.to_list())
