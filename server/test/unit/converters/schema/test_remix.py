import json
import unittest
import unittest.mock

import anndata
import numpy
import pandas as pd

from server.converters.schema import remix


class TestFieldParsing(unittest.TestCase):

    def test_is_curie(self):
        self.assertTrue(remix.is_curie("EFO:00001"))
        self.assertTrue(remix.is_curie("UBERON:123456"))
        self.assertTrue(remix.is_curie("HsapDv:0001"))
        self.assertFalse(remix.is_curie("UBERON"))
        self.assertFalse(remix.is_curie("UBERON:"))
        self.assertFalse(remix.is_curie("123456"))

    def test_is_ontology_field(self):
        self.assertTrue(remix.is_ontology_field("tissue_ontology_term_id"))
        self.assertTrue(remix.is_ontology_field("cell_type_ontology_term_id"))
        self.assertFalse(remix.is_ontology_field("cell_ontology"))
        self.assertFalse(remix.is_ontology_field("method"))

    def test_get_label_field_name(self):
        self.assertEqual("tissue", remix.get_label_field_name("tissue_ontology_term_id"))
        self.assertEqual("cell_type", remix.get_label_field_name("cell_type_ontology_term_id"))

    def test_split_suffix(self):
        self.assertEqual(("UBERON:1234", " (organoid)"), remix.split_suffix("UBERON:1234 (organoid)"))
        self.assertEqual(("UBERON:1234", " (cell culture)"), remix.split_suffix("UBERON:1234 (cell culture)"))
        self.assertEqual(("UBERON:1234", ""), remix.split_suffix("UBERON:1234"))
        self.assertEqual(("UBERON:1234 (something)", ""), remix.split_suffix("UBERON:1234 (something)"))

    @unittest.mock.patch("server.converters.schema.ontology.get_ontology_label")
    def test_get_curie_and_label(self, mock_get_ontology_label):
        mock_get_ontology_label.return_value = "test label"
        self.assertEqual(
            remix.get_curie_and_label("UBERON:1234"),
            ("UBERON:1234", "test label")
        )
        self.assertEqual(
            remix.get_curie_and_label("UBERON:1234 (cell culture)"),
            ("UBERON:1234 (cell culture)", "test label (cell culture)")
        )
        self.assertEqual(
            remix.get_curie_and_label("whatever"),
            ("", "whatever")
        )


class TestManipulateAnndata(unittest.TestCase):

    def setUp(self):

        self.cell_count = 20
        self.gene_count = 200
        X = numpy.random.randint(0, 1000, (self.cell_count, self.gene_count))
        uns = {"organism": "monkey", "experiment": "monkey experiment"}
        obs = pd.DataFrame(
            index=[f"Cell{d}" for d in range(self.cell_count)],
            columns=["tissue", "CellType"],
            data=[["lung", "epithelial"]] * (self.cell_count // 2) + [["lung", "endothelial"]] * (self.cell_count // 2)
        )
        var = pd.DataFrame(index=[f"SEPT{d}" for d in range(self.gene_count)])

        self.adata = anndata.AnnData(X=X, obs=obs, var=var, uns=uns)

    def test_safe_add_field(self):

        remix.safe_add_field(self.adata.obs, "tissue", ["monkey lung"] * self.cell_count)
        self.assertEqual(self.adata.obs["tissue_original"].tolist(), ["lung"] * self.cell_count)
        self.assertEqual(self.adata.obs["tissue"].tolist(), ["monkey lung"] * self.cell_count)

        remix.safe_add_field(self.adata.uns, "contributors", [{"name": "contributor1"}, {"name": "contributor2"}])
        self.assertEqual(
            self.adata.uns["contributors"],
            json.dumps([{"name": "contributor1"}, {"name": "contributor2"}])
        )
