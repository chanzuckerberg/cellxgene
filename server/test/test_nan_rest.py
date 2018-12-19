from http import HTTPStatus
from subprocess import Popen
import unittest
import time

import requests

LOCAL_URL = "http://127.0.0.1:5005/"
VERSION = "v0.2"
URL_BASE = f"{LOCAL_URL}api/{VERSION}/"

BAD_FILTER = {"filter": {"obs": {"annotation_value": [{"name": "xyz"}]}}}


class WithNaNs(unittest.TestCase):
    """Test Case for endpoints"""

    @classmethod
    def setUpClass(cls):
        cls.ps = Popen(
            ["cellxgene", "launch", "server/test/test_datasets/nan.h5ad", "--debug"]
        )
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

    def test_errors(self):
        endpoints = ["annotations/obs", "annotations/var", "data/obs", "data/var"]
        for endpoint in endpoints:
            url = f"{URL_BASE}{endpoint}"
            result = self.session.get(url)
            self.assertEqual(result.status_code, HTTPStatus.INTERNAL_SERVER_ERROR)


class WithoutNaNs(unittest.TestCase):
    """Test Case for endpoints"""

    @classmethod
    def setUpClass(cls):
        cls.ps = Popen(
            [
                "cellxgene",
                "launch",
                "server/test/test_datasets/nan.h5ad",
                "--nan-to-num",
                "--debug",
            ]
        )
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

    def test_errors(self):
        endpoints = ["annotations/obs", "annotations/var", "data/obs", "data/var"]
        for endpoint in endpoints:
            url = f"{URL_BASE}{endpoint}"
            result = self.session.get(url)
            self.assertEqual(result.status_code, HTTPStatus.OK)
