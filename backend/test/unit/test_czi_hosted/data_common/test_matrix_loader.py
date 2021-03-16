import os
import shutil
import tempfile
import time
import unittest

from backend.czi_hosted.common.config.app_config import AppConfig
from backend.common.errors import DatasetAccessError
from backend.czi_hosted.data_common.matrix_loader import MatrixDataCacheManager
from backend.test.unit.test_czi_hosted import FIXTURES_ROOT


class MatrixCacheTest(unittest.TestCase):
    def setup(self):
        pass

    def make_temporay_datasets(self, dirname, num):
        source = f"{FIXTURES_ROOT}/pbmc3k.cxg"
        for i in range(num):
            target = os.path.join(dirname, str(i) + ".cxg")
            shutil.copytree(source, target)

    def use_dataset(self, matrix_cache, dirname, app_config, dataset_index):
        with matrix_cache.data_adaptor(None, os.path.join(dirname, str(dataset_index) + ".cxg"), app_config) as adaptor:
            pass
        return adaptor

    def use_dataset_with_error(self, matrix_cache, dirname, app_config, dataset_index):
        try:
            with matrix_cache.data_adaptor(None, os.path.join(dirname, str(dataset_index) + ".cxg"), app_config):
                raise DatasetAccessError("something bad happened")
        except DatasetAccessError:
            # the MatrixDataCacheManager rethrows the exception, so catch and ignore
            pass

    def get_datasets(self, matrix_cache, dirname):
        datasets = matrix_cache.datasets
        result = {}
        for k, v in datasets.items():
            # filter out the dirname and the .cxg from the name
            newk = int(k[1][len(dirname) + 1 : -4])
            result[newk] = v

        return result

    def check_datasets(self, matrix_cache, dirname, expected):
        res = self.get_datasets(matrix_cache, dirname)
        actual = res.keys()
        self.assertSetEqual(set(actual), set(expected))

    def test_basic(self):
        with tempfile.TemporaryDirectory() as dirname:
            self.make_temporay_datasets(dirname, 5)
            app_config = AppConfig()
            m = MatrixDataCacheManager(max_cached=3, timelimit_s=None)

            # should have only dataset 0
            self.use_dataset(m, dirname, app_config, 0)
            self.check_datasets(m, dirname, [0])

            # should have datasets 0, 1
            self.use_dataset(m, dirname, app_config, 1)
            self.check_datasets(m, dirname, [0, 1])

            # should have datasets 0, 1, 2
            self.use_dataset(m, dirname, app_config, 2)
            self.check_datasets(m, dirname, [0, 1, 2])

            # should have datasets 1, 2, 3
            self.use_dataset(m, dirname, app_config, 3)
            self.check_datasets(m, dirname, [1, 2, 3])

            # use dataset 1, making is more recent than dataset 2
            self.use_dataset(m, dirname, app_config, 1)
            self.check_datasets(m, dirname, [1, 2, 3])

            # use dataset 4, should have 1,3,4
            self.use_dataset(m, dirname, app_config, 4)
            self.check_datasets(m, dirname, [1, 3, 4])

            # use dataset 4 a few more times, get the count to 3
            self.use_dataset(m, dirname, app_config, 4)
            self.use_dataset(m, dirname, app_config, 4)

            datasets = self.get_datasets(m, dirname)
            self.assertEqual(datasets[1].num_access, 2)
            self.assertEqual(datasets[3].num_access, 1)
            self.assertEqual(datasets[4].num_access, 3)

    def test_timelimit(self):
        with tempfile.TemporaryDirectory() as dirname:
            self.make_temporay_datasets(dirname, 2)

            app_config = AppConfig()
            m = MatrixDataCacheManager(max_cached=3, timelimit_s=1)

            adaptor = self.use_dataset(m, dirname, app_config, 0)
            adaptor1 = self.use_dataset(m, dirname, app_config, 0)
            self.assertTrue(adaptor is adaptor1)

            # wait until the timelimit expires and check that there is a new adaptor
            time.sleep(1.1)
            adaptor2 = self.use_dataset(m, dirname, app_config, 0)
            self.assertTrue(adaptor is not adaptor2)
            self.check_datasets(m, dirname, [0])

            # now load a different dataset and see if dataset 0 gets evicted
            time.sleep(1.1)
            self.use_dataset(m, dirname, app_config, 1)
            self.check_datasets(m, dirname, [1])

    def test_access_error(self):
        with tempfile.TemporaryDirectory() as dirname:
            self.make_temporay_datasets(dirname, 1)

            app_config = AppConfig()
            m = MatrixDataCacheManager(max_cached=3, timelimit_s=1)

            # use the 0 datasets
            self.use_dataset(m, dirname, app_config, 0)
            self.check_datasets(m, dirname, [0])

            # use the 0 datasets, but this time a DatasetAccessError is raised.
            # verify that dataset is removed from the cache.
            self.use_dataset_with_error(m, dirname, app_config, 0)
            self.check_datasets(m, dirname, [])
