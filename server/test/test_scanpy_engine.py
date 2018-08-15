import unittest

from server.app.scanpy_engine.scanpy_engine import ScanpyEngine


class UtilTest(unittest.TestCase):
    def setUp(self):
        self.data = ScanpyEngine("example-dataset/", schema="data_schema.json")

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
        filter = {
            "filter": {
                "var": {
                    "index": [1, 99, [200, 300]]
                },
                "obs": {
                    "index": [1, 99, [1000, 2000]]
                }
            }
        }
        data = self.data.filter_dataframe(filter["filter"])
        self.assertEqual(data.shape, (1002, 102))


if __name__ == '__main__':
    unittest.main()
