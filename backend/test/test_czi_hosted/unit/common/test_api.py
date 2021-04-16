import shutil
import time
import unittest
import zlib
from http import HTTPStatus
import hashlib

import pandas as pd
import requests

from backend.czi_hosted.data_common.matrix_loader import MatrixDataType
from backend.test.test_czi_hosted.unit import (
    data_with_tmp_annotations,
    make_fbs,
    start_test_server,
    stop_test_server,
)
from backend.test import PROJECT_ROOT, FIXTURES_ROOT, decode_fbs
from backend.test.fixtures.fixtures import pbmc3k_colors

BAD_FILTER = {"filter": {"obs": {"annotation_value": [{"name": "xyz"}]}}}


# TODO (mweiden): remove ANNOTATIONS_ENABLED and Annotation subclasses when annotations are no longer experimental


class EndPoints(object):
    ANNOTATIONS_ENABLED = True

    def test_initialize(self):
        endpoint = "schema"
        url = f"{self.URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = result.json()
        self.assertEqual(result_data["schema"]["dataframe"]["nObs"], 2638)
        self.assertEqual(len(result_data["schema"]["annotations"]["obs"]), 2)
        self.assertEqual(
            len(result_data["schema"]["annotations"]["obs"]["columns"]), 6 if self.ANNOTATIONS_ENABLED else 5
        )

    def test_config(self):
        endpoint = "config"
        url = f"{self.URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = result.json()
        self.assertIn("library_versions", result_data["config"])
        self.assertEqual(result_data["config"]["displayNames"]["dataset"], "pbmc3k")

    def test_get_layout_fbs(self):
        endpoint = "layout/obs"
        url = f"{self.URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 8)
        self.assertIsNotNone(df["columns"])
        self.assertSetEqual(
            set(df["col_idx"]),
            {"pca_0", "pca_1", "tsne_0", "tsne_1", "umap_0", "umap_1", "draw_graph_fr_0", "draw_graph_fr_1"},
        )
        self.assertIsNone(df["row_idx"])
        self.assertEqual(len(df["columns"]), df["n_cols"])

    def test_put_layout_fbs(self):
        # first check that re-embedding is turned on
        result = self.session.get(f"{self.URL_BASE}config")
        config_data = result.json()
        re_embed = config_data["config"]["parameters"]["enable-reembedding"]
        if not re_embed:
            return
        # attempt to reembed with umap over 100 cells.
        endpoint = "layout/obs"
        url = f"{self.URL_BASE}{endpoint}"
        data = {}
        data["filter"] = {}
        data["filter"]["obs"] = {}
        data["filter"]["obs"]["index"] = list(range(100))
        data["method"] = "umap"
        result = self.session.put(url, json=data)

        self.assertEqual(result.status_code, HTTPStatus.OK)
        result_data = result.json()
        self.assertIsInstance(result_data, dict)
        self.assertEqual(result_data["type"], "float32")
        self.assertTrue(result_data["name"].startswith("reembed:umap_"))
        self.assertIsInstance(result_data["dims"], list)
        self.assertEqual(len(result_data["dims"]), 2)
        dims = result_data["dims"]
        self.assertTrue(dims[0].startswith("reembed:umap_") and dims[0].endswith("_0"))
        self.assertTrue(dims[1].startswith("reembed:umap_") and dims[1].endswith("_1"))

    def test_bad_filter(self):
        endpoint = "data/var"
        url = f"{self.URL_BASE}{endpoint}"
        result = self.session.put(url, json=BAD_FILTER)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

    def test_get_annotations_obs_fbs(self):
        endpoint = "annotations/obs"
        url = f"{self.URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 6 if self.ANNOTATIONS_ENABLED else 5)
        self.assertIsNotNone(df["columns"])
        self.assertIsNone(df["row_idx"])
        self.assertEqual(len(df["columns"]), df["n_cols"])
        obs_index_col_name = self.schema["schema"]["annotations"]["obs"]["index"]
        self.assertCountEqual(
            df["col_idx"],
            [obs_index_col_name, "n_genes", "percent_mito", "n_counts", "louvain"]
            + (["cluster-test"] if self.ANNOTATIONS_ENABLED else []),
        )

    def test_get_annotations_obs_keys_fbs(self):
        endpoint = "annotations/obs"
        query = "annotation-name=n_genes&annotation-name=percent_mito"
        url = f"{self.URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 2)
        self.assertIsNotNone(df["columns"])
        self.assertIsNone(df["row_idx"])
        self.assertEqual(len(df["columns"]), df["n_cols"])
        self.assertCountEqual(df["col_idx"], ["n_genes", "percent_mito"])

    def test_get_annotations_obs_error(self):
        endpoint = "annotations/obs"
        query = "annotation-name=notakey"
        url = f"{self.URL_BASE}{endpoint}?{query}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

    def test_diff_exp(self):
        endpoint = "diffexp/obs"
        url = f"{self.URL_BASE}{endpoint}"
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
        url = f"{self.URL_BASE}{endpoint}"
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
        url = f"{self.URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 1838)
        self.assertEqual(df["n_cols"], 2)
        self.assertIsNotNone(df["columns"])
        self.assertIsNone(df["row_idx"])
        self.assertEqual(len(df["columns"]), df["n_cols"])
        var_index_col_name = self.schema["schema"]["annotations"]["var"]["index"]
        self.assertCountEqual(df["col_idx"], [var_index_col_name, "n_cells"])

    def test_get_annotations_var_keys_fbs(self):
        endpoint = "annotations/var"
        query = "annotation-name=n_cells"
        url = f"{self.URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 1838)
        self.assertEqual(df["n_cols"], 1)
        self.assertIsNotNone(df["columns"])
        self.assertIsNone(df["row_idx"])
        self.assertEqual(len(df["columns"]), df["n_cols"])
        self.assertCountEqual(df["col_idx"], ["n_cells"])

    def test_get_annotations_var_error(self):
        endpoint = "annotations/var"
        query = "annotation-name=notakey"
        url = f"{self.URL_BASE}{endpoint}?{query}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

    def test_data_mimetype_error(self):
        endpoint = "data/var"
        header = {"Accept": "xxx"}
        url = f"{self.URL_BASE}{endpoint}"
        result = self.session.put(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.NOT_ACCEPTABLE)

    def test_fbs_default(self):
        endpoint = "data/var"
        url = f"{self.URL_BASE}{endpoint}"
        result = self.session.put(url)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

        filter = {"filter": {"var": {"index": [0, 1, 4]}}}
        result = self.session.put(url, json=filter)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")

    def test_data_put_fbs(self):
        endpoint = "data/var"
        url = f"{self.URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.put(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

    def test_data_get_fbs(self):
        endpoint = "data/var"
        url = f"{self.URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

    def test_data_put_filter_fbs(self):
        endpoint = "data/var"
        url = f"{self.URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        filter = {"filter": {"var": {"index": [0, 1, 4]}}}
        result = self.session.put(url, headers=header, json=filter)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 3)
        self.assertIsNotNone(df["columns"])
        self.assertIsNone(df["row_idx"])
        self.assertEqual(len(df["columns"]), df["n_cols"])
        self.assertListEqual(df["col_idx"].tolist(), [0, 1, 4])

    def test_data_get_filter_fbs(self):
        index_col_name = self.schema["schema"]["annotations"]["var"]["index"]
        endpoint = "data/var"
        query = f"var:{index_col_name}=SIK1"
        url = f"{self.URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 1)

    def test_data_get_unknown_filter_fbs(self):
        index_col_name = self.schema["schema"]["annotations"]["var"]["index"]
        endpoint = "data/var"
        query = f"var:{index_col_name}=UNKNOWN"
        url = f"{self.URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 0)

    def test_data_put_single_var(self):
        endpoint = "data/var"
        url = f"{self.URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        index_col_name = self.schema["schema"]["annotations"]["var"]["index"]
        var_filter = {"filter": {"var": {"annotation_value": [{"name": index_col_name, "values": ["RER1"]}]}}}
        result = self.session.put(url, headers=header, json=var_filter)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 1)

    def test_colors(self):
        endpoint = "colors"
        url = f"{self.URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = result.json()
        self.assertEqual(result_data, pbmc3k_colors)

    def test_static(self):
        endpoint = "static"
        file = "assets/favicon.ico"
        url = f"{self.server}/{endpoint}/{file}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)

    def test_genesets_config(self):
        result = self.session.get(f"{self.URL_BASE}config")
        config_data = result.json()
        params = config_data["config"]["parameters"]
        annotations_genesets = params["annotations_genesets"]
        annotations_genesets_readonly = params["annotations_genesets_readonly"]
        annotations_genesets_summary_methods = params["annotations_genesets_summary_methods"]
        self.assertTrue(annotations_genesets)
        self.assertTrue(annotations_genesets_readonly)
        self.assertEqual(annotations_genesets_summary_methods, ["mean"])

    def test_get_genesets(self):
        endpoint = "genesets"
        url = f"{self.URL_BASE}{endpoint}"
        result = self.session.get(url, headers={"Accept": "application/json"})
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = result.json()
        self.assertIsNotNone(result_data["genesets"])

    def test_get_summaryvar(self):
        index_col_name = self.schema["schema"]["annotations"]["var"]["index"]
        endpoint = "summarize/var"

        # single column
        filter = f"var:{index_col_name}=F5"
        query = f"method=mean&{filter}"
        query_hash = hashlib.sha1(query.encode()).hexdigest()
        url = f"{self.URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 1)
        self.assertEqual(df["col_idx"], [query_hash])
        self.assertAlmostEqual(df["columns"][0][0], -0.110451095)

        # multi-column
        col_names = ["F5", "BEB3", "SIK1"]
        filter = "&".join([f"var:{index_col_name}={name}" for name in col_names])
        query = f"method=mean&{filter}"
        query_hash = hashlib.sha1(query.encode()).hexdigest()
        url = f"{self.URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 1)
        self.assertEqual(df["col_idx"], [query_hash])
        self.assertAlmostEqual(df["columns"][0][0], -0.16628358)

    def test_post_summaryvar(self):
        index_col_name = self.schema["schema"]["annotations"]["var"]["index"]
        endpoint = "summarize/var"
        headers = {"Content-Type": "application/x-www-form-urlencoded", "Accept": "application/octet-stream"}

        # single column
        filter = f"var:{index_col_name}=F5"
        query = f"method=mean&{filter}"
        query_hash = hashlib.sha1(query.encode()).hexdigest()
        url = f"{self.URL_BASE}{endpoint}?key={query_hash}"
        result = self.session.post(url, headers=headers, data=query)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 1)
        self.assertEqual(df["col_idx"], [query_hash])
        self.assertAlmostEqual(df["columns"][0][0], -0.110451095)

        # multi-column
        col_names = ["F5", "BEB3", "SIK1"]
        filter = "&".join([f"var:{index_col_name}={name}" for name in col_names])
        query = f"method=mean&{filter}"
        query_hash = hashlib.sha1(query.encode()).hexdigest()
        url = f"{self.URL_BASE}{endpoint}?key={query_hash}"
        result = self.session.post(url, headers=headers, data=query)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 1)
        self.assertEqual(df["col_idx"], [query_hash])
        self.assertAlmostEqual(df["columns"][0][0], -0.16628358)

    def _setupClass(child_class, command_line):
        child_class.ps, child_class.server = start_test_server(command_line)
        child_class.URL_BASE = f"{child_class.server}/api/v0.2/"
        child_class.session = requests.Session()
        for i in range(90):
            try:
                result = child_class.session.get(f"{child_class.URL_BASE}schema")
                child_class.schema = result.json()
            except requests.exceptions.ConnectionError:
                time.sleep(1)


class EndPointsAnnotations(EndPoints):
    def test_get_schema_existing_writable(self):
        self._test_get_schema_writable("cluster-test")

    def test_get_user_annotations_existing_obs_keys_fbs(self):
        self._test_get_user_annotations_obs_keys_fbs(
            "cluster-test",
            {"unassigned", "one", "two", "three", "four", "five", "six", "seven"},
        )

    def test_put_user_annotations_obs_fbs(self):
        endpoint = "annotations/obs"
        query = "annotation-collection-name=test_annotations"
        url = f"{self.URL_BASE}{endpoint}?{query}"
        n_rows = self.data.get_shape()[0]
        fbs = make_fbs({"cat_A": pd.Series(["label_A"] * n_rows, dtype="category")})
        result = self.session.put(url, data=zlib.compress(fbs))
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        self.assertEqual(result.json(), {"status": "OK"})
        self._test_get_schema_writable("cat_A")
        self._test_get_user_annotations_obs_keys_fbs("cat_A", {"label_A"})

    def _test_get_user_annotations_obs_keys_fbs(self, annotation_name, columns):
        endpoint = "annotations/obs"
        query = f"annotation-name={annotation_name}"
        url = f"{self.URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.content)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 1)
        self.assertListEqual(df["col_idx"], [annotation_name])
        self.assertEqual(set(df["columns"][0]), columns)
        self.assertIsNone(df["row_idx"])
        self.assertEqual(len(df["columns"]), df["n_cols"])

    def _test_get_schema_writable(self, cluster_name):
        endpoint = "schema"
        url = f"{self.URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = result.json()
        columns = result_data["schema"]["annotations"]["obs"]["columns"]
        matching_columns = [c for c in columns if c["name"] == cluster_name]
        self.assertEqual(len(matching_columns), 1)
        self.assertTrue(matching_columns[0]["writable"])


class EndPointsAnndata(unittest.TestCase, EndPoints):
    """Test Case for endpoints"""

    ANNOTATIONS_ENABLED = False

    @classmethod
    def setUpClass(cls):
        cls._setupClass(
            cls,
            [
                f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad",
                "--disable-annotations",
                "--experimental-enable-reembedding",
            ],
        )

    @classmethod
    def tearDownClass(cls):
        stop_test_server(cls.ps)

    @property
    def annotations_enabled(self):
        return False


class EndPointsCxg(unittest.TestCase, EndPoints):
    """Test Case for endpoints"""

    ANNOTATIONS_ENABLED = False

    @classmethod
    def setUpClass(cls):
        cls._setupClass(cls, [f"{FIXTURES_ROOT}/pbmc3k.cxg", "--disable-annotations"])

    @classmethod
    def tearDownClass(cls):
        stop_test_server(cls.ps)


class EndPointsAnndataAnnotations(unittest.TestCase, EndPointsAnnotations):
    """Test Case for endpoints"""

    ANNOTATIONS_ENABLED = True

    @classmethod
    def setUpClass(cls):
        cls.data, cls.tmp_dir, cls.annotations = data_with_tmp_annotations(
            MatrixDataType.H5AD, annotations_fixture=True
        )
        cls._setupClass(cls, ["--annotations-file", cls.annotations.output_file, cls.data.get_location()])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp_dir)
        stop_test_server(cls.ps)


class EndPointsCxgAnnotations(unittest.TestCase, EndPointsAnnotations):
    """Test Case for endpoints"""

    ANNOTATIONS_ENABLED = True

    @classmethod
    def setUpClass(cls):
        cls.data, cls.tmp_dir, cls.annotations = data_with_tmp_annotations(MatrixDataType.CXG, annotations_fixture=True)
        cls._setupClass(cls, ["--annotations-file", cls.annotations.output_file, cls.data.get_location()])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp_dir)
        stop_test_server(cls.ps)
