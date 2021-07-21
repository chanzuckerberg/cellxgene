import unittest

import numpy as np

from backend.common.compute import diffexp_generic
from backend.server.data_common.matrix_loader import MatrixDataLoader
from backend.test.test_server.unit import app_config
from backend.test import PROJECT_ROOT


class DiffExpTest(unittest.TestCase):
    """Tests the diffexp returns the expected results for one test case, using the h5ad
    adaptor types and different algorithms."""

    def load_dataset(self, path, extra_server_config={}, extra_dataset_config={}):
        config = app_config(path, extra_server_config=extra_server_config, extra_dataset_config=extra_dataset_config)
        loader = MatrixDataLoader(path)
        adaptor = loader.open(config)
        return adaptor

    def get_mask(self, adaptor, start, stride):
        """Simple function to return a mask or rows"""
        rows = adaptor.get_shape()[0]
        sel = list(range(start, rows, stride))
        mask = np.zeros(rows, dtype=bool)
        mask[sel] = True
        return mask

    def compare_diffexp_results(self, results, expects):
        self.assertEqual(len(results), len(expects))
        for result, expect in zip(results, expects):
            self.assertEqual(result[0], expect[0])
            self.assertTrue(np.isclose(result[1], expect[1], 1e-6, 1e-4))
            self.assertTrue(np.isclose(result[2], expect[2], 1e-6, 1e-4))
            self.assertTrue(np.isclose(result[3], expect[3], 1e-6, 1e-4))

    def check_1_10_2_10(self, results):
        """Checks the results for a specific set of rows selections"""

        positive_expects = [
            [1712, -0.5525154, 0.0051788902660723345, 1.0],
            [1575, 1.0317602, 0.007830310753043345, 1.0],
            [693, 0.4703904, 0.008715846769131548, 1.0],
            [916, 0.9567287, 0.009080596532247588, 1.0],
            [77, 0.02665649, 0.010070392939027756, 1.0],
            [782, -1.0981874, 0.010161745218916036, 1.0],
            [913, 0.5683986, 0.010782030711612685, 1.0],
            [910, 0.83164597, 0.014596411069229197, 1.0],
            [1727, 0.4127781, 0.015168372104237176, 1.0],
            [1443, -0.8241895, 0.015337080567465522, 1.0],
        ]
        negative_expects = [
            [956, 0.016060986, 0.0008649321884808977, 1.0],
            [1124, 0.96602094, 0.0011717216548271284, 1.0],
            [1809, 1.1110606, 0.0019304405196777848, 1.0],
            [1754, 0.5201581, 0.005691734062127954, 1.0],
            [948, 1.6390722, 0.006622111055981219, 1.0],
            [1810, 0.78618884, 0.007055917428377063, 1.0],
            [779, 1.5241305, 0.007202934422407284, 1.0],
            [576, 0.97873515, 0.008272092578813124, 1.0],
            [538, 0.89114505, 0.01062259019889307, 1.0],
            [436, 0.3119122, 0.01127515110543434, 1.0],
        ]

        self.compare_diffexp_results(results["positive"], positive_expects)
        self.compare_diffexp_results(results["negative"], negative_expects)

    def get_X_col(self, adaptor, cols):
        varmask = np.zeros(adaptor.get_shape()[1], dtype=bool)
        varmask[cols] = True
        return adaptor.get_X_array(None, varmask)

    def test_anndata_default(self):
        """Test an anndata adaptor with its default diffexp algorithm (diffexp_generic)"""
        adaptor = self.load_dataset(f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad")
        maskA = self.get_mask(adaptor, 1, 10)
        maskB = self.get_mask(adaptor, 2, 10)
        results = adaptor.compute_diffexp_ttest(maskA, maskB, 10)
        self.check_1_10_2_10(results)


def test_h5ad_default(self):
    """Test a h5ad adaptor with its default diffexp algorithm (diffexp_cxg)"""
    adaptor = self.load_dataset(f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad")
    maskA = self.get_mask(adaptor, 1, 10)
    maskB = self.get_mask(adaptor, 2, 10)

    # run it through the adaptor
    results = adaptor.compute_diffexp_ttest(maskA, maskB, 10)
    self.check_1_10_2_10(results)

    # run it directly

    results = diffexp_generic.diffexp_ttest(adaptor, maskA, maskB, 10)
    self.check_1_10_2_10(results)


def test_h5ad_generic(self):
    """Test a h5ad adaptor with the generic adaptor"""
    adaptor = self.load_dataset(f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad")
    maskA = self.get_mask(adaptor, 1, 10)
    maskB = self.get_mask(adaptor, 2, 10)
    # run it directly
    results = diffexp_generic.diffexp_ttest(adaptor, maskA, maskB, 10)
    self.check_1_10_2_10(results)
