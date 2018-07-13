import unittest

from server.app.scanpy_engine.scanpy_engine import ScanpyEngine


class UtilTest(unittest.TestCase):
    def setUp(self):
        self.data = ScanpyEngine("example-dataset/", schema="data_schema.json")

    def test_init(self):
        assert self.data.cell_count == 2638
        assert self.data.gene_count == 1838
        epsilon = 0.000005
        assert self.data.data.X[0,0] - -0.17146951 < epsilon

    def test_schema(self):
        assert self.data.schema == {'CellName': {'type': 'string', 'variabletype': 'categorical', 'displayname': 'Name', 'include': True}, 'n_genes': {'type': 'int', 'variabletype': 'continuous', 'displayname': 'Num Genes', 'include': True}, 'percent_mito': {'type': 'float', 'variabletype': 'continuous', 'displayname': 'Mitochondrial Percentage', 'include': True}, 'n_counts': {'type': 'float', 'variabletype': 'continuous', 'displayname': 'Num Counts', 'include': True}, 'louvain': {'type': 'string', 'variabletype': 'categorical', 'displayname': 'Louvain Cluster', 'include': True}}

    def test_cells(self):
        cells = self.data.cells()
        assert "AAACATACAACCAC-1" in cells
        assert len(cells) == 2638

    def test_genes(self):
        genes = self.data.genes()
        assert "SEPT4" in genes
        assert len(genes) == 1838

    def test_filter_categorical(self):
        filter = {"louvain": {"variable_type": "categorical", "value_type": "string", "query": ["B cells"]}}
        filtered_data = self.data.filter_cells(filter)
        assert filtered_data.shape == (342, 1838)
        louvain_vals = filtered_data.obs['louvain'].tolist()
        assert "B cells" in louvain_vals
        assert "NK cells" not in louvain_vals

    def test_filter_continuous(self):
        # print(self.data.data.obs["n_genes"].tolist())
        filter = {"n_genes": {"variable_type": "continuous", "value_type": "int", "query": {"min": 300, "max": 400}}}
        filtered_data = self.data.filter_cells(filter)
        assert filtered_data.shape == (71, 1838)
        n_genes_vals = filtered_data.obs['n_genes'].tolist()
        for val in n_genes_vals:
            assert 300 <= val <= 400

    def test_metadata(self):
        metadata = self.data.metadata(df=self.data.data)
        assert len(metadata) == 2638
        assert 'louvain' in metadata[0]

    def test_create_graph(self):
        graph = self.data.create_graph(df=self.data.data)
        assert graph[0][1] == 0.5545382653143183
        assert graph[0][2] == 0.6021833809031731

    def test_diffexp(self):
        diffexp = self.data.diffexp(["AAACATACAACCAC-1", "AACCGATGGTCATG-1"], ["CCGATAGACCTAAG-1", "GGTGGAGAAGTAGA-1"], 0.5, 7)
        assert diffexp["celllist1"]["topgenes"] == ['EBNA1BP2', 'DIAPH1', 'SLC25A11', 'SNRNP27', 'COMMD8', 'COTL1', 'GTF3A']

    def test_expression(self):
        expression = self.data.expression(cells=["AAACATACAACCAC-1"])
        data_exp = self.data.data[["AAACATACAACCAC-1"], :].X
        for idx in range(len(expression["cells"][0]["e"])):
            assert expression["cells"][0]["e"][idx] == data_exp[idx]


if __name__ == '__main__':
    unittest.main()