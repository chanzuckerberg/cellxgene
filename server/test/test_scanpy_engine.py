import json
from os import path
import pytest
import time
import unittest

import numpy as np
from pandas import Series


from server.app.scanpy_engine.scanpy_engine import ScanpyEngine


class UtilTest(unittest.TestCase):
    def setUp(self):
        self.data = ScanpyEngine("example-dataset/")
        self.data._create_schema()

    def test_init(self):
        self.assertEqual(self.data.cell_count, 2638)
        self.assertEqual(self.data.gene_count, 1838)
        epsilon = 0.000005
        self.assertTrue(self.data.data.X[0, 0] - -0.17146951 < epsilon)

    def test_schema(self):
        self.assertEqual(self.data.schema, {'CellName': {'type': 'string', 'variabletype': 'categorical', 'displayname': 'Name', 'include': True}, 'n_genes': {'type': 'int', 'variabletype': 'continuous', 'displayname': 'Num Genes', 'include': True}, 'percent_mito': {'type': 'float', 'variabletype': 'continuous', 'displayname': 'Mitochondrial Percentage', 'include': True}, 'n_counts': {'type': 'float', 'variabletype': 'continuous', 'displayname': 'Num Counts', 'include': True}, 'louvain': {'type': 'string', 'variabletype': 'categorical', 'displayname': 'Louvain Cluster', 'include': True}})

    def test_cells(self):
        cells = self.data.cells()
        self.assertIn("AAACATACAACCAC-1", cells)
        self.assertEqual(len(cells), 2638)

    def test_genes(self):
        genes = self.data.genes()
        self.assertIn("SEPT4", genes)
        self.assertEqual(len(genes), 1838)

    def test_filter_categorical(self):
        filter = {"louvain": {"variable_type": "categorical", "value_type": "string", "query": ["B cells"]}}
        filtered_data = self.data.filter_cells(filter)
        self.assertEqual(filtered_data.shape, (342, 1838))
        louvain_vals = filtered_data.obs['louvain'].tolist()
        self.assertIn("B cells", louvain_vals)
        self.assertNotIn("NK cells", louvain_vals)

    def test_filter_continuous(self):
        # print(self.data.data.obs["n_genes"].tolist())
        filter = {"n_genes": {"variable_type": "continuous", "value_type": "int", "query": {"min": 300, "max": 400}}}
        filtered_data = self.data.filter_cells(filter)
        self.assertEqual(filtered_data.shape, (71, 1838))
        n_genes_vals = filtered_data.obs['n_genes'].tolist()
        for val in n_genes_vals:
            self.assertTrue(300 <= val <= 400)

    def test_metadata(self):
        metadata = self.data.metadata(df=self.data.data)
        self.assertEqual(len(metadata), 2638)
        self.assertIn('louvain', metadata[0])

    @unittest.skip("Umap not producing the same graph on different systems, even with the same seed. Skipping for now")
    def test_create_graph(self):
        graph = self.data.create_graph(df=self.data.data)
        self.assertEqual(graph[0][1], 0.5545382653143183)
        self.assertEqual(graph[0][2], 0.6021833809031731)

    def test_diffexp(self):
        diffexp = self.data.diffexp(["AAACATACAACCAC-1", "AACCGATGGTCATG-1"], ["CCGATAGACCTAAG-1", "GGTGGAGAAGTAGA-1"], 0.5, 7)
        self.assertEqual(diffexp["celllist1"]["topgenes"], ['EBNA1BP2', 'DIAPH1', 'SLC25A11', 'SNRNP27', 'COMMD8', 'COTL1', 'GTF3A'])

    def test_expression(self):
        expression = self.data.expression(cells=["AAACATACAACCAC-1"])
        data_exp = self.data.data[["AAACATACAACCAC-1"], :].X
        for idx in range(len(expression["cells"][0]["e"])):
            self.assertEqual(expression["cells"][0]["e"][idx], data_exp[idx])

    def test_filter_idx(self):
        filter_ = {
            "filter": {
                "var": {
                    "index": [1, 99, [200, 300]]
                },
                "obs": {
                    "index": [1, 99, [1000, 2000]]
                }
            }
        }
        data = self.data.filter_dataframe(filter_["filter"])
        self.assertEqual(data.shape, (1002, 102))

    def test_filter_annotation(self):
        filter_ = {
            "filter": {
                "obs": {
                    "annotation_value": [
                        {"name": "louvain", "values": ["NK cells", "CD8 T cells"]},
                    ]
                }
            }
        }
        data = self.data.filter_dataframe(filter_["filter"])
        self.assertEqual(data.shape, (470, 1838))
        filter_ = {
            "filter": {
                "obs": {
                    "annotation_value": [
                        {"name": "n_counts", "min": 3000},
                    ]
                }
            }
        }
        data = self.data.filter_dataframe(filter_["filter"])
        self.assertEqual(data.shape, (497, 1838))

    def test_filter_complex(self):
        filter_ = {
            "filter": {
                "var": {
                    "index": [1, 99, [200, 300]]
                },
                "obs": {
                    "annotation_value": [
                        {"name": "louvain", "values": ["NK cells", "CD8 T cells"]},
                        {"name": "n_counts", "min": 3000},
                    ],
                    "index": [1, 99, [1000, 2000]]
                }
            }
        }
        data = self.data.filter_dataframe(filter_["filter"])
        self.assertEqual(data.shape, (15, 102))

    def test_obs_and_var_names(self):
        self.assertEqual(np.sum(self.data.data.var["name"].isna()), 0)
        self.assertEqual(np.sum(self.data.data.obs["name"].isna()), 0)

    def test_schema(self):
        with open(path.join(path.dirname(__file__), "schema.json")) as fh:
            schema = json.load(fh)
            self.assertEqual(self.data.schema, schema)

    def test_schema_produces_error(self):
        self.data.data.obs["time"] = Series(list([time.time() for i in range(self.data.cell_count)]),
                                            dtype="datetime64[ns]")
        with pytest.raises(TypeError):
            self.data._create_schema()

    if __name__ == '__main__':
        unittest.main()
