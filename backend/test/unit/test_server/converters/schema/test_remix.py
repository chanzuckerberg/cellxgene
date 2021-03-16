import json
import os
import unittest
import unittest.mock

import anndata
import numpy
import pandas as pd
import scanpy as sc

from backend.server.converters.schema import remix
from backend.test.unit.test_server import PROJECT_ROOT


class TestApplySchema(unittest.TestCase):

    def setUp(self):
        self.source_h5ad_path = f"{PROJECT_ROOT}/backend/test/fixtures/pbmc3k-CSC-gz.h5ad"
        self.output_h5ad_path = f"{PROJECT_ROOT}/backend/test/fixtures/test_remix.h5ad"
        self.config_path = f"{PROJECT_ROOT}/backend/test/fixtures/test_config.yaml"
        self.bad_config_path = f"{PROJECT_ROOT}/backend/test/fixtures/test_bad_config.yaml"

    def tearDown(self):
        try:
            os.remove(self.output_h5ad_path)
        except OSError:
            pass

    @unittest.mock.patch("backend.server.converters.schema.ontology.get_ontology_label")
    def test_apply_schema(self, mock_get_ontology_label):
        mock_get_ontology_label.return_value = "test label"
        remix.apply_schema(self.source_h5ad_path, self.config_path, self.output_h5ad_path)
        new_adata = sc.read_h5ad(self.output_h5ad_path)

        self.assertIn("cell_type", new_adata.obs.columns)
        self.assertListEqual(["test label"], new_adata.obs["cell_type"].unique().tolist())
        self.assertListEqual(
            ["CL:00001", "CL:00002", "CL:00003", "CL:00004", "CL:00005", "CL:00006", "CL:00007", "CL:00008"],
            sorted(new_adata.obs["cell_type_ontology_term_id"].unique().tolist())
        )

        self.assertIn("version", new_adata.uns_keys())

    @unittest.mock.patch("backend.server.converters.schema.ontology.get_ontology_label")
    def test_apply_bad_schema(self, mock_get_ontology_label):
        mock_get_ontology_label.return_value = "test label"
        remix.apply_schema(self.source_h5ad_path, self.bad_config_path, self.output_h5ad_path)
        new_adata = sc.read_h5ad(self.output_h5ad_path)

        # Should refuse to write the version
        self.assertNotIn("version", new_adata.uns_keys())

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

    @unittest.mock.patch("backend.server.converters.schema.ontology.get_ontology_label")
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

    @unittest.mock.patch("backend.server.converters.schema.ontology.get_ontology_label")
    def test_remix_uns(self, mock_get_ontology_label):
        mock_get_ontology_label.return_value = "Pan troglodytes"
        uns_config = {
            "version": {
                "corpora_schema_version": "1.0.0",
                "corpora_encoding_version": "0.1.0"
            },
            "organism_ontology_term_id": "NCBITaxon:9598",
            "contributors": [
                {
                    "name": "scientist",
                    "email": "scientist@science.com"
                }
            ]
        }

        remix.remix_uns(self.adata, uns_config)

        self.assertEqual(
            sorted(self.adata.uns_keys()),
            sorted(["organism_original", "organism", "organism_ontology_term_id",
                    "contributors", "version", "experiment"])
        )

        self.assertEqual(self.adata.uns['organism'], "Pan troglodytes")
        self.assertEqual(self.adata.uns['organism_original'], "monkey")
        self.assertEqual(self.adata.uns['organism_ontology_term_id'], "NCBITaxon:9598")
        self.assertEqual(self.adata.uns['contributors'],
                         json.dumps([{"name": "scientist", "email": "scientist@science.com"}]))

    @unittest.mock.patch("backend.server.converters.schema.ontology.get_ontology_label")
    def test_remix_obs(self, mock_get_ontology_label):
        mock_get_ontology_label.return_value = "lung (in a monkey)"
        obs_config = {
            "tissue_ontology_term_id": {
                "tissue": {
                    "lung": "UBERON:00000"
                }
            },
            "cell_color": {
                "CellType": {
                    "epithelial": "fuschia",
                    "endothelial": "khaki"
                }
            },
            "sex": "male"
        }

        remix.remix_obs(self.adata, obs_config)
        self.assertEqual(
            sorted(self.adata.obs_keys()),
            sorted(["tissue", "tissue_ontology_term_id", "tissue_original", "CellType", "cell_color", "sex"])
        )

        self.assertTrue(all(v == "lung" for v in self.adata.obs.tissue_original))
        self.assertTrue(all(v == "UBERON:00000" for v in self.adata.obs.tissue_ontology_term_id))
        self.assertTrue(all(v == "lung (in a monkey)" for v in self.adata.obs.tissue))
        self.assertTrue(all(v == "male" for v in self.adata.obs.sex))
        self.assertTrue(all(v in (("epithelial", "fuschia"), ("endothelial", "khaki"))
                            for v in zip(self.adata.obs.CellType, self.adata.obs.cell_color)))


class TestFixupGeneSymbols(unittest.TestCase):

    def setUp(self):
        self.seurat_path = f"{PROJECT_ROOT}/server/test/fixtures/schema_test_data/seurat_tutorial.h5ad"
        self.seurat_merged_path = f"{PROJECT_ROOT}/server/test/fixtures/schema_test_data/seurat_tutorial_merged.h5ad"
        self.sctransform_path = f"{PROJECT_ROOT}/server/test/fixtures/schema_test_data/sctransform.h5ad"
        self.sctransform_merged_path = f"{PROJECT_ROOT}/server/test/fixtures/schema_test_data/sctransform_merged.h5ad"

        # There's lots of MALAT1, but it doesn't collide with any other names,
        # so it shouldn't change during merging.
        self.stable_gene = "MALAT1"

    def test_fixup_gene_symbols_seurat(self):

        if not os.path.isfile(self.seurat_path):
            return unittest.skip(
                "Skipping gene symbol conversion tests because test h5ads are not present. To create them, "
                "run server/test/fixtures/schema_test_data/generate_test_data.sh"
            )

        original_adata = sc.read_h5ad(self.seurat_path)
        merged_adata = sc.read_h5ad(self.seurat_merged_path)

        fixup_config = {"X": "log1p", "counts": "raw", "scale.data": "log1p"}

        fixed_adata = remix.fixup_gene_symbols(original_adata, fixup_config)

        self.assertEqual(
            merged_adata.layers["counts"][:, merged_adata.var.index == self.stable_gene].sum(),
            fixed_adata.raw.X[:, fixed_adata.var.index == self.stable_gene].sum()
        )
        self.assertAlmostEqual(
            merged_adata.X[:, merged_adata.var.index == self.stable_gene].sum(),
            fixed_adata.X[:, fixed_adata.var.index == self.stable_gene].sum()
        )

        self.assertAlmostEqual(
            merged_adata.layers["scale.data"][:, merged_adata.var.index == self.stable_gene].sum(),
            fixed_adata.layers["scale.data"][:, fixed_adata.var.index == self.stable_gene].sum()
        )

    def test_fixup_gene_symbols_sctransform(self):

        if not os.path.isfile(self.sctransform_path):
            return unittest.skip(
                "Skipping gene symbol conversion tests because test h5ads are not present. To create them, "
                "run server/test/fixtures/schema_test_data/generate_test_data.sh"
            )

        original_adata = sc.read_h5ad(self.sctransform_path)
        merged_adata = sc.read_h5ad(self.sctransform_merged_path)

        fixup_config = {"X": "log1p", "counts": "raw"}

        fixed_adata = remix.fixup_gene_symbols(original_adata, fixup_config)

        # sctransform does a bunch of stuff, including slightly modifying the
        # raw counts. So we can't assert for exact equality the way we do with
        # the vanilla seurat tutorial. But, the results should still be very
        # close.
        merged_raw_stable = merged_adata.layers["counts"][:, merged_adata.var.index == self.stable_gene].sum()
        fixed_raw_stable = fixed_adata.raw.X[:, fixed_adata.var.index == self.stable_gene].sum()
        self.assertLess(abs(merged_raw_stable - fixed_raw_stable), .001 * merged_raw_stable)

        self.assertAlmostEqual(
            merged_adata.X[:, merged_adata.var.index == self.stable_gene].sum(),
            fixed_adata.X[:, fixed_adata.var.index == self.stable_gene].sum(),
            0
        )
