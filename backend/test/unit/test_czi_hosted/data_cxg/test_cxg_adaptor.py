import unittest

from backend.common.utils.data_locator import DataLocator
from backend.czi_hosted.data_cxg.cxg_adaptor import CxgAdaptor
from backend.test.unit.test_czi_hosted import app_config
from backend.test.unit import FIXTURES_ROOT
from backend.test.fixtures.fixtures import pbmc3k_colors


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
