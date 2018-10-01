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

    def test_get_annotations_obs(self):
        endpoint = "annotations/obs"
        url = f"{URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["names"], ["n_genes", "percent_mito", "n_counts", "louvain", "name"])
        self.assertEqual(len(result_data["data"]), 2638)
        self.assertEqual(len(result_data["data"][0]), 6)

    def test_get_annotations_obs_keys(self):
        endpoint = "annotations/obs"
        query = "annotation-name=n_genes&annotation-name=percent_mito"
        url = f"{URL_BASE}{endpoint}?{query}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["names"], ["n_genes", "percent_mito"])
        self.assertEqual(len(result_data["data"][0]), 3)

    def test_get_annotations_obs_error(self):
        endpoint = "annotations/obs"
        query = "annotation-name=notakey"
        url = f"{URL_BASE}{endpoint}?{query}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 404)

    def test_put_annotations_obs(self):
        endpoint = "annotations/obs"
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
        self.assertEqual(result_data["names"], ["n_genes", "percent_mito", "n_counts", "louvain", "name"])
        self.assertEqual(len(result_data["data"]), 15)

    def test_filter_put_annotations_obs(self):
        endpoint = "annotations/obs"
        query = "annotation-name=n_genes&annotation-name=percent_mito"
        url = f"{URL_BASE}{endpoint}?{query}"
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
        self.assertEqual(result_data["names"], ["n_genes", "percent_mito"])
        self.assertEqual(len(result_data["data"][0]), 3)
        self.assertEqual(len(result_data["data"]), 15)

    def test_diff_exp(self):
        endpoint = "diffexp/obs"
        url = f"{URL_BASE}{endpoint}"
        params = {
            "mode": "topN",
            "set1": {
                "filter": {
                    "obs": {"annotation_value": [
                        {"name": "louvain", "values": ["NK cells"]}
                    ]
                    }
                }
            },
            "set2": {
                "filter": {
                    "obs": {"annotation_value": [
                        {"name": "louvain", "values": ["CD8 T cells"]}
                    ]
                    }
                }
            },
            "count": 7
        }
        result = self.session.post(url, json=params)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(len(result_data), 7)

    def test_diff_exp_indices(self):
        endpoint = "diffexp/obs"
        url = f"{URL_BASE}{endpoint}"
        params = {
            "mode": "topN",
            "set1": {
                "filter": {
                    "obs": {
                        "index": [[0, 500]]
                    }
                }
            },
            "set2": {
                "filter": {
                    "obs": {
                        "index": [[500, 1000]]
                    }
                }
            }
        }
        result = self.session.post(url, json=params)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(len(result_data), 10)

    def test_get_annotations_var(self):
        endpoint = "annotations/var"
        url = f"{URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["names"], ["n_cells", "name"])
        self.assertEqual(len(result_data["data"]), 1838)
        self.assertEqual(len(result_data["data"][0]), 3)

    def test_get_annotations_var_keys(self):
        endpoint = "annotations/var"
        query = "annotation-name=n_cells"
        url = f"{URL_BASE}{endpoint}?{query}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["names"], ["n_cells"])
        self.assertEqual(len(result_data["data"][0]), 2)

    def test_get_annotations_var_error(self):
        endpoint = "annotations/var"
        query = "annotation-name=notakey"
        url = f"{URL_BASE}{endpoint}?{query}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 404)

    def test_put_annotations_var(self):
        endpoint = "annotations/var"
        url = f"{URL_BASE}{endpoint}"
        var_filter = {
            "filter": {
                "var": {
                    "annotation_value": [
                        {"name": "name", "values": ["ATAD3C", "RER1"]},
                    ]
                }
            }
        }
        result = self.session.put(url, json=var_filter)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["names"], ["n_cells", "name"])
        self.assertEqual(len(result_data["data"]), 2)

    def test_filter_put_annotations_var(self):
        endpoint = "annotations/var"
        query = "annotation-name=n_cells"
        url = f"{URL_BASE}{endpoint}?{query}"
        var_filter = {
            "filter": {
                "var": {
                    "annotation_value": [
                        {"name": "name", "values": ["ATAD3C", "RER1"]},
                    ]
                }
            }
        }
        result = self.session.put(url, json=var_filter)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["names"], ["n_cells"])
        self.assertEqual(len(result_data["data"][0]), 2)
        self.assertEqual(len(result_data["data"]), 2)

    def test_get_data(self):
        endpoint = "data/obs"
        query = "accept-type=application/json"
        url = f"{URL_BASE}{endpoint}?{query}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(len(result_data["obs"]), 2638)

    def test_data_mimetype_error(self):
        endpoint = "data/obs"
        query = "accept-type=xxx"
        url = f"{URL_BASE}{endpoint}?{query}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 406)
        # no accept type
        url = f"{URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 406)

    def test_data_filter(self):
        endpoint = "data/obs"
        query = "accept-type=application/json&obs:louvain=NK cells&obs:louvain=CD8 T cells&obs:n_counts=3000,*"
        url = f"{URL_BASE}{endpoint}?{query}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(len(result_data["obs"]), 38)

    def test_data_put(self):
        endpoint = "data/obs"
        url = f"{URL_BASE}{endpoint}"
        header = {"Accept": "application/json"}
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
        result = self.session.put(url, headers=header, json=obs_filter)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(len(result_data["obs"]), 15)

    def test_data_put_single_var(self):
        endpoint = "data/obs"
        url = f"{URL_BASE}{endpoint}"
        header = {"Accept": "application/json"}
        var_filter = {
            "filter": {
                "var": {
                    "annotation_value": [
                        {"name": "name", "values": ["RER1"]},
                    ]
                }
            }
        }
        result = self.session.put(url, headers=header, json=var_filter)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(len(result_data["obs"][0]), 2)

    def test_static(self):
        endpoint = "static"
        file = "js/service-worker.js"
        url = f"{LOCAL_URL}{endpoint}/{file}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
