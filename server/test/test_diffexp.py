import unittest
import subprocess
import sys
import shutil
from server.data_common.matrix_loader import MatrixDataLoader
from server.common.app_config import AppConfig
import server.compute.diffexp_cxg as diffexp_cxg
import server.compute.diffexp_generic as diffexp_generic
from server.compute.diffexp_generic import compute_indices_hash_version, compute_indices_hash
import numpy as np

from server.test import PROJECT_ROOT


class DiffExpTest(unittest.TestCase):
    """Tests the diffexp returns the expected results for one test case, using different
    adaptor types and different algorithms."""

    def setUp(self):
        diffexp = f"{PROJECT_ROOT}/server/test/test_datasets/pbmc3k.cxg/diffexp"
        shutil.rmtree(diffexp, ignore_errors=True)

    def tearnDown(self):
        diffexp = f"{PROJECT_ROOT}/server/test/test_datasets/pbmc3k.cxg/diffexp"
        shutil.rmtree(diffexp, ignore_errors=True)

    def load_dataset(self, path):
        app_config = AppConfig()
        app_config.single_dataset__datapath = path
        app_config.server__verbose = True
        app_config.complete_config()
        loader = MatrixDataLoader(path)
        adaptor = loader.open(app_config)
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
            self.assertAlmostEqual(result[1], expect[1])
            self.assertAlmostEqual(result[2], expect[2])
            self.assertAlmostEqual(result[3], expect[3])

    def check_1_10_2_10(self, results):
        """Checks the results for a specific set of rows selections"""
        expects = [
            [956, 0.016060986, 0.0008649321884808977, 1.0],
            [1124, 0.96602094, 0.0011717216548271284, 1.0],
            [1809, 1.1110606, 0.0019304405196777848, 1.0],
            [1712, -0.5525154, 0.0051788902660723345, 1.0],
            [1754, 0.5201581, 0.005691734062127954, 1.0],
            [948, 1.6390722, 0.006622111055981219, 1.0],
            [1810, 0.78618884, 0.007055917428377063, 1.0],
            [779, 1.5241305, 0.007202934422407284, 1.0],
            [1575, 1.0317602, 0.007830310753043345, 1.0],
            [576, 0.97873515, 0.008272092578813124, 1.0],
        ]
        self.compare_diffexp_results(results, expects)

    def test_anndata_default(self):
        """Test an anndata adaptor with its default diffexp algorithm (diffexp_generic)"""
        adaptor = self.load_dataset(f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad")
        maskA = self.get_mask(adaptor, 1, 10)
        maskB = self.get_mask(adaptor, 2, 10)
        results = adaptor.compute_diffexp_ttest(maskA, maskB, 10)
        self.check_1_10_2_10(results)

    def test_cxg_default(self):
        """Test a cxg adaptor with its default diffexp algorithm (diffexp_cxg)"""
        adaptor = self.load_dataset(f"{PROJECT_ROOT}/server/test/test_datasets/pbmc3k.cxg")
        maskA = self.get_mask(adaptor, 1, 10)
        maskB = self.get_mask(adaptor, 2, 10)

        # run it through the adaptor
        results = adaptor.compute_diffexp_ttest(maskA, maskB, 10)
        self.check_1_10_2_10(results)

        # run it directly
        results = diffexp_cxg.diffexp_ttest(adaptor, maskA, maskB, 10)
        self.check_1_10_2_10(results)

    def test_cxg_generic(self):
        """Test a cxg adaptor with the generic adaptor"""
        adaptor = self.load_dataset(f"{PROJECT_ROOT}/server/test/test_datasets/pbmc3k.cxg")
        maskA = self.get_mask(adaptor, 1, 10)
        maskB = self.get_mask(adaptor, 2, 10)
        # run it directly
        results = diffexp_generic.diffexp_ttest(adaptor, maskA, maskB, 10)
        self.check_1_10_2_10(results)

    def test_cache(self):

        dataset = f"{PROJECT_ROOT}/server/test/test_datasets/pbmc3k.cxg"
        script = f"{PROJECT_ROOT}/server/converters/diffexp_precompute.py"
        subprocess.check_call(f"{sys.executable} {script} {dataset} --overwrite", shell=True)
        adaptor = self.load_dataset(dataset)
        obsvals = adaptor.query_obs_array("louvain")[:]
        filterA = np.where(obsvals == "CD8 T cells")[0]
        filterB = np.where(obsvals == "CD4 T cells")[0]
        rows = adaptor.get_shape()[0]
        maskA = np.zeros(rows, dtype=bool)
        maskA[filterA] = True
        maskB = np.zeros(rows, dtype=bool)
        maskB[filterB] = True

        res_generic = diffexp_generic.diffexp_ttest(adaptor, maskA, maskB, 10)

        # run it directly
        cache_hash_version = compute_indices_hash_version()
        res_cache = diffexp_cxg.diffexp_ttest(adaptor, maskA, maskB, 10, cache_hash_version=cache_hash_version)
        self.compare_diffexp_results(res_cache, res_generic)
        res_no_cache = diffexp_cxg.diffexp_ttest(adaptor, maskA, maskB, 10, cache_hash_version=None)
        self.compare_diffexp_results(res_no_cache, res_generic)

        # run it through the adaptor
        res_cache_adaptor = adaptor.compute_diffexp_ttest(maskA, maskB, 10, use_cache=True)
        self.compare_diffexp_results(res_cache_adaptor, res_generic)
        res_no_cache_adaptor = adaptor.compute_diffexp_ttest(maskA, maskB, 10, use_cache=False)
        self.compare_diffexp_results(res_no_cache_adaptor, res_generic)

        # check that the cache works
        hashA = compute_indices_hash(filterA, cache_hash_version)
        assert adaptor.diffexp_cache(hashA) == 3
        assert diffexp_cxg.query_cache(adaptor, filterA, cache_hash_version) is not None
        hashB = compute_indices_hash(filterB, cache_hash_version)
        assert adaptor.diffexp_cache(hashB) == 2
        assert diffexp_cxg.query_cache(adaptor, filterB, cache_hash_version) is not None
        hashN = compute_indices_hash(filterA[:-1], cache_hash_version)
        assert adaptor.diffexp_cache(hashN) is None
        assert diffexp_cxg.query_cache(adaptor, filterA[:-1], cache_hash_version) is None

        # check that a new cache hash version does not fail
        res_cache = diffexp_cxg.diffexp_ttest(adaptor, maskA, maskB, 10, cache_hash_version=cache_hash_version + 1)
        self.compare_diffexp_results(res_cache, res_generic)
