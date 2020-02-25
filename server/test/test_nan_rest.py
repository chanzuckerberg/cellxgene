from http import HTTPStatus
from subprocess import Popen
import unittest
import time
import math

import server.test.decode_fbs as decode_fbs

import requests

LOCAL_URL = "http://127.0.0.1:5006/"
VERSION = "v0.2"
URL_BASE = f"{LOCAL_URL}api/{VERSION}/"

BAD_FILTER = {"filter": {"obs": {"annotation_value": [{"name": "xyz"}]}}}


class WithNaNs(unittest.TestCase):
    """Test Case for endpoints"""

    @classmethod
    def setUpClass(cls):
        cls.ps = Popen(["cellxgene", "launch", "test/test_datasets/nan.h5ad", "--verbose", "--port", "5006"])
        session = requests.Session()
        for i in range(90):
            try:
                session.get(f"{URL_BASE}schema")
            except requests.exceptions.ConnectionError:
                time.sleep(1)

    @classmethod
    def tearDownClass(cls):
        try:
            cls.ps.terminate()
        except ProcessLookupError:
            pass

    def setUp(self):
        self.session = requests.Session()

    def test_initialize(self):
        endpoint = "schema"
        url = f"{URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)

    def test_data(self):
        endpoint = "data/var"
        url = f"{URL_BASE}{endpoint}"
        result = self.session.put(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertTrue(math.isnan(df["columns"][3][3]))

    def test_annotation_obs(self):
        endpoint = "annotations/obs"
        url = f"{URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertTrue(math.isnan(df["columns"][2][0]))

    def test_annotation_var(self):
        endpoint = "annotations/var"
        url = f"{URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertTrue(math.isnan(df["columns"][2][0]))
