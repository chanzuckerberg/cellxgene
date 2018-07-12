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



if __name__ == '__main__':
    unittest.main()