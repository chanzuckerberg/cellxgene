from http import HTTPStatus
from subprocess import Popen
import unittest
import time

import requests

import decode_fbs

LOCAL_URL = "http://127.0.0.1:5005/"
VERSION = "v0.2"
URL_BASE = f"{LOCAL_URL}api/{VERSION}/"

BAD_FILTER = {"filter": {"obs": {"annotation_value": [{"name": "xyz"}]}}}


class EndPoints(unittest.TestCase):
    """Test Case for endpoints"""

    @classmethod
    def setUpClass(cls):
        cls.ps = Popen(["cellxgene", "launch", "example-dataset/pbmc3k.h5ad", "--debug", "--port", "5005"])
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
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = result.json()
        self.assertEqual(result_data["schema"]["dataframe"]["nObs"], 2638)
        self.assertEqual(len(result_data["schema"]["annotations"]["obs"]), 5)

    def test_config(self):
        endpoint = "config"
        url = f"{URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = result.json()
        self.assertIn("library_versions", result_data["config"])
        self.assertEqual(result_data["config"]["displayNames"]["dataset"], "pbmc3k")
        self.assertEqual(len(result_data["config"]["features"]), 4)

    def test_get_layout_fbs(self):
        endpoint = "layout/obs"
        url = f"{URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df['n_rows'], 2638)
        self.assertEqual(df['n_cols'], 2)
        self.assertIsNotNone(df['columns'])
        self.assertIsNone(df['col_idx'])
        self.assertIsNone(df['row_idx'])
        self.assertEqual(len(df['columns']), df['n_cols'])

    def test_bad_filter(self):
        endpoint = "data/var"
        url = f"{URL_BASE}{endpoint}"
        result = self.session.put(url, json=BAD_FILTER)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

    def test_get_annotations_obs_fbs(self):
        endpoint = "annotations/obs"
        url = f"{URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df['n_rows'], 2638)
        self.assertEqual(df['n_cols'], 5)
        self.assertIsNotNone(df['columns'])
        self.assertIsNotNone(df['col_idx'])
        self.assertIsNone(df['row_idx'])
        self.assertEqual(len(df['columns']), df['n_cols'])
        self.assertListEqual(df['col_idx'], ['name', 'n_genes', 'percent_mito', 'n_counts', 'louvain'])

    def test_get_annotations_obs_keys_fbs(self):
        endpoint = "annotations/obs"
        query = "annotation-name=n_genes&annotation-name=percent_mito"
        url = f"{URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df['n_rows'], 2638)
        self.assertEqual(df['n_cols'], 2)
        self.assertIsNotNone(df['columns'])
        self.assertIsNotNone(df['col_idx'])
        self.assertIsNone(df['row_idx'])
        self.assertEqual(len(df['columns']), df['n_cols'])
        self.assertListEqual(df['col_idx'], ['n_genes', 'percent_mito'])

    def test_get_annotations_obs_error(self):
        endpoint = "annotations/obs"
        query = "annotation-name=notakey"
        url = f"{URL_BASE}{endpoint}?{query}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

    def test_diff_exp(self):
        endpoint = "diffexp/obs"
        url = f"{URL_BASE}{endpoint}"
        params = {
            "mode": "topN",
            "set1": {"filter": {"obs": {"annotation_value": [{"name": "louvain", "values": ["NK cells"]}]}}},
            "set2": {"filter": {"obs": {"annotation_value": [{"name": "louvain", "values": ["CD8 T cells"]}]}}},
            "count": 7,
        }
        result = self.session.post(url, json=params)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = result.json()
        self.assertEqual(len(result_data), 7)

    def test_diff_exp_indices(self):
        endpoint = "diffexp/obs"
        url = f"{URL_BASE}{endpoint}"
        params = {
            "mode": "topN",
            "count": 10,
            "set1": {"filter": {"obs": {"index": [[0, 500]]}}},
            "set2": {"filter": {"obs": {"index": [[500, 1000]]}}},
        }
        result = self.session.post(url, json=params)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = result.json()
        self.assertEqual(len(result_data), 10)

    def test_get_annotations_var_fbs(self):
        endpoint = "annotations/var"
        url = f"{URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df['n_rows'], 1838)
        self.assertEqual(df['n_cols'], 2)
        self.assertIsNotNone(df['columns'])
        self.assertIsNotNone(df['col_idx'])
        self.assertIsNone(df['row_idx'])
        self.assertEqual(len(df['columns']), df['n_cols'])
        self.assertListEqual(df['col_idx'], ['name', 'n_cells'])

    def test_get_annotations_var_keys_fbs(self):
        endpoint = "annotations/var"
        query = "annotation-name=n_cells"
        url = f"{URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df['n_rows'], 1838)
        self.assertEqual(df['n_cols'], 1)
        self.assertIsNotNone(df['columns'])
        self.assertIsNotNone(df['col_idx'])
        self.assertIsNone(df['row_idx'])
        self.assertEqual(len(df['columns']), df['n_cols'])
        self.assertListEqual(df['col_idx'], ['n_cells'])

    def test_get_annotations_var_error(self):
        endpoint = "annotations/var"
        query = "annotation-name=notakey"
        url = f"{URL_BASE}{endpoint}?{query}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

    def test_data_mimetype_error(self):
        endpoint = f"data/var"
        header = {"Accept": "xxx"}
        url = f"{URL_BASE}{endpoint}"
        result = self.session.put(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.NOT_ACCEPTABLE)

    def test_fbs_default(self):
        endpoint = f"data/var"
        url = f"{URL_BASE}{endpoint}"
        result = self.session.put(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")

    def test_data_put_fbs(self):
        endpoint = f"data/var"
        url = f"{URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.put(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df['n_rows'], 2638)
        self.assertEqual(df['n_cols'], 1838)
        self.assertIsNotNone(df['columns'])
        self.assertListEqual(df['col_idx'].tolist(), [])
        self.assertIsNone(df['row_idx'])
        self.assertEqual(len(df['columns']), df['n_cols'])

    def test_data_put_filter_fbs(self):
        endpoint = f"data/var"
        url = f"{URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        filter = {
            "filter": {
                "var": {
                    "index": [0, 1, 4]
                }
            }
        }
        result = self.session.put(url, headers=header, json=filter)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df['n_rows'], 2638)
        self.assertEqual(df['n_cols'], 3)
        self.assertIsNotNone(df['columns'])
        self.assertIsNotNone(df['col_idx'])
        self.assertIsNone(df['row_idx'])
        self.assertEqual(len(df['columns']), df['n_cols'])
        self.assertListEqual(df['col_idx'].tolist(), [0, 1, 4])

    def test_data_put_single_var(self):
        endpoint = f"data/var"
        url = f"{URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        var_filter = {"filter": {"var": {"annotation_value": [{"name": "name", "values": ["RER1"]}]}}}
        result = self.session.put(url, headers=header, json=var_filter)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 1)

    def test_static(self):
        endpoint = "static"
        file = "js/service-worker.js"
        url = f"{LOCAL_URL}{endpoint}/{file}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
