from http import HTTPStatus
import math

import backend.test.decode_fbs as decode_fbs

from backend.czi_hosted.common.config.app_config import AppConfig
from backend.test import FIXTURES_ROOT
from backend.test.test_czi_hosted.unit import BaseTest

VERSION = "v0.2"
BAD_FILTER = {"filter": {"obs": {"annotation_value": [{"name": "xyz"}]}}}


class WithNaNs(BaseTest):
    """Test Case for endpoints"""

    @classmethod
    def setUpClass(cls):
        app_config = AppConfig()
        app_config.update_server_config(single_dataset__datapath=f"{FIXTURES_ROOT}/nan.h5ad")
        app_config.update_default_dataset_config(user_annotations__enable=True)
        super().setUpClass(app_config)
        cls.app.testing = True
        cls.client = cls.app.test_client()

    def setUp(self):
        self.session = self.client
        self.url_base = "api/v0.2/"

    def test_initialize(self):
        endpoint = "schema"
        url = f"{self.url_base}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)

    def test_data(self):
        endpoint = "data/var"
        url = f"{self.url_base}{endpoint}"
        filter = {"filter": {"var": {"index": [[0, 20]]}}}
        header = {"Accept": "application/octet-stream"}
        result = self.session.put(url, headers=header, json=filter)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertTrue(math.isnan(df["columns"][3][3]))

    def test_annotation_obs(self):
        endpoint = "annotations/obs"
        url = f"{self.url_base}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertTrue(math.isnan(df["columns"][2][0]))

    def test_annotation_var(self):
        endpoint = "annotations/var"
        url = f"{self.url_base}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertTrue(math.isnan(df["columns"][2][0]))
