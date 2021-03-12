import unittest

from backend.server.common.data_locator import DataLocator
from backend.server.data_cxg.cxg_adaptor import CxgAdaptor
from backend.server.test import FIXTURES_ROOT, app_config
from backend.server.test.fixtures.fixtures import pbmc3k_colors


class TestCxgAdaptor(unittest.TestCase):
    def test_get_colors(self):
        data = self.get_data("pbmc3k.cxg")
        self.assertDictEqual(data.get_colors(), pbmc3k_colors)
        data = self.get_data("pbmc3k_v0.cxg")
        self.assertDictEqual(data.get_colors(), dict())

    def get_data(self, fixture):
        data_locator = f"{FIXTURES_ROOT}/{fixture}"
        config = app_config(data_locator)
        return CxgAdaptor(DataLocator(data_locator), config)
