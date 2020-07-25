import shutil
import unittest
import json

import anndata

from server.common.data_locator import DataLocator
from server.converters.cxgtool import write_cxg, create_cxg_group_metadata
from server.data_cxg.cxg_adaptor import CxgAdaptor
from server.test import PROJECT_ROOT, app_config, random_string
from server.test.test_datasets.fixtures import pbmc3k_colors


class TestCxgAdaptor(unittest.TestCase):
    def setUp(self) -> None:
        self.fixtures = []

    def tearDown(self) -> None:
        try:
            for data_locator in self.fixtures:
                print("REMOVING ", data_locator)
                shutil.rmtree(data_locator)
        except FileNotFoundError:
            pass

    def test_cxg_category_colors(self):
        data = self.convert_pbmc3k(extract_colors=True)
        self.assertEqual(data.get_colors(), pbmc3k_colors)
        data = self.convert_pbmc3k(extract_colors=False)
        self.assertEqual(data.get_colors(), {})

    def convert_pbmc3k(self, **kwargs):
        rand_str = random_string(8)
        data_locator = f"/tmp/test_{rand_str}.cxg"
        self.fixtures.append(data_locator)
        source_h5ad = anndata.read_h5ad(f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad")
        cxg_group_metadata = create_cxg_group_metadata(
            adata=source_h5ad,
            basefname="pbmc3k.h5ad",
            title="pbmc3k",
            **kwargs
        )
        write_cxg(
            adata=source_h5ad,
            container=data_locator,
            cxg_group_metadata=cxg_group_metadata
        )
        config = app_config(data_locator)
        return CxgAdaptor(DataLocator(data_locator), config)
