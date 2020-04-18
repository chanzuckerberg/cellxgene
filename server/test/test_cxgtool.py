import random
import shutil
import string
import unittest

import anndata

from server.common.data_locator import DataLocator
from server.converters.cxgtool import write_cxg
from server.data_cxg.cxg_adaptor import CxgAdaptor
from server.test import PROJECT_ROOT, app_config
from server.test.test_datasets.fixtures import pbmc3k_colors


class TestCxgAdaptor(unittest.TestCase):

    def setUp(self) -> None:
        random_string = ''.join(random.choice(string.ascii_letters) for _ in range(8))
        self.data_locator = f"/tmp/test_{random_string}.cxg"
        self.source_h5ad = anndata.read_h5ad(f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad")
        write_cxg(
            adata=self.source_h5ad,
            container=self.data_locator,
            title="pbmc3k",
            extract_colors=True
        )
        config = app_config(self.data_locator)
        self.data = CxgAdaptor(DataLocator(self.data_locator), config)

    def tearDown(self) -> None:
        try:
            shutil.rmtree(self.data_locator)
        except FileNotFoundError:
            pass

    def test_cxg_category_colors(self):
        self.assertEqual(self.data.get_colors(), pbmc3k_colors)
