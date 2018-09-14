import requests
from subprocess import Popen
import unittest
import time

LOCAL_URL = "http://127.0.0.1:5005/"
VERSION = "v0.2"
URL_BASE = f"{LOCAL_URL}api/{VERSION}/"



class EndPoints(unittest.TestCase):
    """Test Case for endpoints"""

    @classmethod
    def setUpClass(cls):
        cls.ps = Popen(["cellxgene", "scanpy", "example-dataset/"])
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
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["schema"]["dataframe"]["nObs"], 2638)
        self.assertEqual(len(result_data["schema"]["annotations"]["obs"]), 5)

    def test_config(self):
        endpoint = "config"
        url = f"{URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["config"]["displayNames"]["dataset"], "example-dataset")
        self.assertEqual(len(result_data["config"]["features"]), 4)

    def test_get_layout(self):
        endpoint = "layout/obs"
        url = f"{URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["layout"]["ndims"], 2)
        self.assertEqual(len(result_data["layout"]["coordinates"]), 2638)

    def test_put_layout(self):
        endpoint = "layout/obs"
        url = f"{URL_BASE}{endpoint}"
        obs_filter = {
            "filter": {
                "obs": {
                    "annotation_value": [
                        {"name": "louvain", "values": ["NK cells", "CD8 T cells"]},
                        {"name": "n_counts", "min": 3000},
                    ],
                    "index": [1, 99, [1000, 2000]]
                }
            }
        }
        result = self.session.put(url, json=obs_filter)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(len(result_data["layout"]["coordinates"]), 15)

    def test_get_annotations(self):
        endpoint = "annotations/obs"
        url = f"{URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["names"], ["n_genes", "percent_mito", "n_counts", "louvain", "name"])
        self.assertEqual(len(result_data["data"]), 2638)
        self.assertEqual(len(result_data["data"][0]), 6)

    def test_get_annotations_keys(self):
        endpoint = "annotations/obs"
        query = "annotation-name=n_genes&annotation-name=percent_mito"
        url = f"{URL_BASE}{endpoint}?{query}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["names"], ["n_genes", "percent_mito"])
        self.assertEqual(len(result_data["data"][0]), 3)

    def test_get_annotations_error(self):
        endpoint = "annotations/obs"
        query = "annotation-name=notakey"
        url = f"{URL_BASE}{endpoint}?{query}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 404)

    def test_static(self):
        endpoint = "static"
        file = "js/service-worker.js"
        url = f"{LOCAL_URL}{endpoint}/{file}"
        result = self.session.get(url)
        print(result.status_code)
        self.assertEqual(result.status_code, 200)
