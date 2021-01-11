import unittest

from server.data_cxg.cxg_adaptor import CxgAdaptor
from server.czi.czi_cxg_adaptor import CZI_CxgAdaptor
from server.common.config.app_config import AppConfig
from server.app.app import get_data_adaptor
from server.test import FIXTURES_ROOT


class TestCziCxgAdaptor(unittest.TestCase):

    def test_config_adator(self):
        """Test the configuration for the CZI_cxg_adaptor"""

        app_config = AppConfig()
        app_config.update_server_config(
            app__flask_secret_key="secret",
            multi_dataset__dataroot=dict(
                cxg=dict(dataroot=FIXTURES_ROOT, base_url="cxg"),
                czi_cxg=dict(dataroot=FIXTURES_ROOT, base_url="czi_cxg"),
            )
        )
        app_config.add_dataroot_config("czi_cxg")
        app_config.update_single_config_from_path_and_value(
            ["per_dataset_config", "czi_cxg", "adaptor", "cxg", "type"], "czi_cxg_adaptor"
        )

        app_config.complete_config()

        # The cxg dataroot should use the default CxgAdaptor
        with get_data_adaptor(app_config, "cxg", "pbmc3k.cxg") as adaptor:
            self.assertEqual(type(adaptor), CxgAdaptor)

        # The czi_cxg dataroot should use the default CZI_CxgAdaptor
        with get_data_adaptor(app_config, "czi_cxg", "pbmc3k.cxg") as adaptor:
            self.assertEqual(type(adaptor), CZI_CxgAdaptor)
