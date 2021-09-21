from http import HTTPStatus
import unittest
import math
from test.unit import start_test_server, stop_test_server
from test import FIXTURES_ROOT
import test.decode_fbs as decode_fbs

import requests

VERSION = "v0.2"
BAD_FILTER = {"filter": {"obs": {"annotation_value": [{"name": "xyz"}]}}}


class WithNaNs(unittest.TestCase):
    """Test Case for endpoints"""

    @classmethod
    def setUpClass(cls):
        cls.ps, cls.server = start_test_server([f"{FIXTURES_ROOT}/nan.h5ad"])

    @classmethod
    def tearDownClass(cls):
        stop_test_server(cls.ps)

    def setUp(self):
        self.session = requests.Session()
        self.url_base = f"{self.server}/api/{VERSION}/"

    def test_initialize(self):
        endpoint = "schema"
        url = f"{self.url_base}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)

    def test_data(self):
        endpoint = "data/var"
        url = f"{self.url_base}{endpoint}"
        filter = {"filter": {"var": {"index": [[0, 20]]}}}
        result = self.session.put(url, json=filter)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertTrue(math.isnan(df["columns"][3][3]))

    def test_annotation_obs(self):
        endpoint = "annotations/obs"
        url = f"{self.url_base}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertTrue(math.isnan(df["columns"][2][0]))

    def test_annotation_var(self):
        endpoint = "annotations/var"
        url = f"{self.url_base}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertTrue(math.isnan(df["columns"][2][0]))
