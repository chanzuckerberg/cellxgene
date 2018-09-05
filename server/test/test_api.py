import unittest
import requests


class EndPoints(unittest.TestCase):
    """Test Case for endpoints"""

    def setUp(self):
        # Local
        self.local_url = "http://127.0.0.1:5005/"
        self.version = "v0.2"
        self.url_base = "{local_url}api/{version}/".format(local_url=self.local_url, version=self.version)
        self.session = requests.Session()

    def test_initialize(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="schema")
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["schema"]["dataframe"]["nObs"], 2638)
        self.assertEqual(len(result_data["schema"]["annotations"]["obs"]), 5)

    def test_config(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="config")
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["config"]["displayNames"]["dataset"], "example-dataset")
        self.assertEqual(len(result_data["config"]["features"]), 4)

    def test_get_layout(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="layout/obs")
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["layout"]["ndims"], 2)
        self.assertEqual(len(result_data["layout"]["coordinates"]), 2638)

    def test_get_annotations(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="annotation/obs")
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["names"], ["n_genes", "percent_mito", "n_counts", "louvain", "name"])
        self.assertEqual(len(result_data["data"]), 2638)
        self.assertEqual(len(result_data["data"][0]), 6)

    def test_get_annotations_keys(self):
        url = "{base}{endpoint}?{query}".format(base=self.url_base, endpoint="annotation/obs",
                                                query="annotation-keys=n_genes,percent_mito")
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
        result_data = result.json()
        self.assertEqual(result_data["names"], ["n_genes", "percent_mito"])
        self.assertEqual(len(result_data["data"][0]), 3)

    def test_get_annotations_error(self):
        url = "{base}{endpoint}?{query}".format(base=self.url_base, endpoint="annotation/obs",
                                                query="annotation-keys=notakey")
        result = self.session.get(url)
        self.assertEqual(result.status_code, 404)

    def test_static(self):
        url = "{url}{endpoint}/{file}".format(url=self.local_url, endpoint="static", file="js/service-worker.js")
        result = self.session.get(url)
        self.assertEqual(result.status_code, 200)
