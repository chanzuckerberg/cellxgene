import os
import tempfile
import unittest
import unicodedata

import anndata
import numpy as np
import pandas as pd

from server.common.utils.data_locator import DataLocator
from server.data_anndata.anndata_adaptor import AnndataAdaptor
from server.common.config.app_config import AppConfig


class UnicodeMetadataTest(unittest.TestCase):
    """
    Test that h5ad files containing non-ASCII Unicode characters in metadata
    (e.g. Greek letters like phi) can be loaded without error.

    Regression test for: https://github.com/chanzuckerberg/cellxgene/issues/2754
    """

    def _create_h5ad_with_unicode_obs(self, obs_labels):
        """Helper to create a minimal h5ad file with the given obs labels."""
        n_obs = len(obs_labels)
        n_var = 5
        X = np.random.rand(n_obs, n_var).astype(np.float32)
        obs = pd.DataFrame(
            {"cell_type": pd.Categorical(obs_labels)},
            index=[f"cell_{i}" for i in range(n_obs)],
        )
        var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_var)])
        # Create a minimal embedding so validation passes
        obsm = {"X_umap": np.random.rand(n_obs, 2).astype(np.float32)}
        adata = anndata.AnnData(X=X, obs=obs, var=var, obsm=obsm)
        return adata

    def _load_h5ad(self, filepath):
        """Load an h5ad file via AnndataAdaptor with a minimal config."""
        locator = DataLocator(filepath)
        config = AppConfig()
        config.update_server_config(
            single_dataset__obs_names=None,
            single_dataset__var_names=None,
            single_dataset__datapath=locator.path,
        )
        config.update_server_config(app__flask_secret_key="secret")
        config.update_dataset_config(
            embeddings__names=["umap"],
            presentation__max_categories=100,
            diffexp__lfc_cutoff=0.01,
        )
        config.complete_config()
        return AnndataAdaptor(locator, config)

    def test_load_h5ad_with_greek_phi_in_metadata(self):
        """Loading h5ad with Greek letter phi in cell labels should succeed."""
        labels = [
            "Interstitial M\u03c6 perivascular",
            "T cell",
            "B cell",
            "Monocyte",
        ]
        adata = self._create_h5ad_with_unicode_obs(labels)
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
            adata.write_h5ad(f.name)
            tmppath = f.name
        try:
            data = self._load_h5ad(tmppath)
            self.assertIsNotNone(data)
            self.assertEqual(data.cell_count, 4)
            self.assertEqual(data.gene_count, 5)
        finally:
            os.unlink(tmppath)

    def test_load_h5ad_with_various_unicode_in_metadata(self):
        """Loading h5ad with various Unicode characters should succeed."""
        labels = [
            "M\u03c6 perivascular",      # Greek phi
            "\u03b1\u03b2 T cell",        # Greek alpha beta
            "Na\u00efve B cell",          # Latin i with diaeresis
            "Th\u00fcringen sample",      # u-umlaut
            "caf\u00e9 sample",           # e-acute
        ]
        adata = self._create_h5ad_with_unicode_obs(labels)
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
            adata.write_h5ad(f.name)
            tmppath = f.name
        try:
            data = self._load_h5ad(tmppath)
            self.assertIsNotNone(data)
            self.assertEqual(data.cell_count, 5)
        finally:
            os.unlink(tmppath)

    def test_unicode_normalization_applied(self):
        """Verify that Unicode strings in obs are NFC-normalized after loading."""
        # NFD form of e-acute (two code points: 'e' + combining acute accent)
        nfd_label = "caf" + "\u0065\u0301" + " sample"
        # NFC form (single code point)
        nfc_label = "caf\u00e9 sample"

        labels = [nfd_label, "T cell", "B cell"]
        adata = self._create_h5ad_with_unicode_obs(labels)
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
            adata.write_h5ad(f.name)
            tmppath = f.name
        try:
            data = self._load_h5ad(tmppath)
            cell_types = data.data.obs["cell_type"].tolist()
            # After normalization, the NFD form should be converted to NFC
            self.assertIn(nfc_label, cell_types)
            for val in cell_types:
                if isinstance(val, str):
                    self.assertEqual(val, unicodedata.normalize("NFC", val))
        finally:
            os.unlink(tmppath)


if __name__ == "__main__":
    unittest.main()
