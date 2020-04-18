import unittest

from server.common.data_locator import DataLocator
from server.data_cxg.cxg_adaptor import CxgAdaptor
from server.test import PROJECT_ROOT, app_config
from server.test.test_datasets.fixtures import pbmc3k_colors


class TestCxgAdaptor(unittest.TestCase):

    def setUp(self):
        data_locator = f"{PROJECT_ROOT}/server/test/test_datasets/pbmc3k.cxg"
        config = app_config(data_locator)
        self.data = CxgAdaptor(DataLocator(data_locator), config)

    def test_get_colors(self):
        self.assertEqual(self.data.get_colors(), pbmc3k_colors)
