import shutil
import time
import unittest
import zlib
from http import HTTPStatus
import tempfile
from os import path
import hashlib
from os.path import basename, splitext

import pandas as pd
import requests
import numpy as np

from parameterized import parameterized_class

import test.decode_fbs as decode_fbs
from server.data_common.matrix_loader import MatrixDataType
from test.unit import (
    data_with_tmp_annotations,
    make_fbs,
    start_test_server,
    stop_test_server,
)
from test.fixtures.fixtures import pbmc3k_colors
from test import PROJECT_ROOT, FIXTURES_ROOT

BAD_FILTER = {"filter": {"obs": {"annotation_value": [{"name": "xyz"}]}}}


# TODO (mweiden): remove ANNOTATIONS_ENABLED and Annotation subclasses when annotations are no longer experimental


class EndPoints(object):
    ANNOTATIONS_ENABLED = True
    GENESETS_READONLY = False

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

        # Check that all schema types are legal
        legal_types = ["boolean", "string", "categorical", "float32", "int32"]
        self.assertEqual(result_data["schema"]["dataframe"]["type"], "float32")
        for column in result_data["schema"]["annotations"]["obs"]["columns"]:
            self.assertIn(column["type"], legal_types)
        for column in result_data["schema"]["annotations"]["var"]["columns"]:
            self.assertIn(column["type"], legal_types)

    def test_config(self):
        endpoint = "config"
        url = f"{self.URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = result.json()
        self.assertIn("library_versions", result_data["config"])

        if hasattr(self, "data_locator"):
            title = splitext(basename(self.data_locator))[0]
        else:
            title = "pbmc3k"
        self.assertEqual(result_data["config"]["displayNames"]["dataset"], title)
        self.assertIsNotNone(result_data["config"]["parameters"])

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
        for column in df["columns"]:
            self.assertEqual(column.dtype, np.float32)

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
        for column in df["columns"]:
            if type(column) is np.ndarray:
                self.assertIn(column.dtype, [np.float32, np.int32])

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
        for column in df["columns"]:
            if type(column) is np.ndarray:
                self.assertIn(column.dtype, [np.float32, np.int32])

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
        for column in df["columns"]:
            if type(column) is np.ndarray:
                self.assertIn(column.dtype, [np.float32, np.int32])

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
        for column in df["columns"]:
            if type(column) is np.ndarray:
                self.assertIn(column.dtype, [np.float32, np.int32])

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
        for column in df["columns"]:
            if type(column) is np.ndarray:
                self.assertIn(column.dtype, [np.float32, np.int32])

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
        self.assertEqual(annotations_genesets_readonly, self.GENESETS_READONLY)
        self.assertEqual(annotations_genesets_summary_methods, ["mean"])

    def test_get_genesets(self):
        endpoint = "genesets"
        url = f"{self.URL_BASE}{endpoint}"
        result = self.session.get(url, headers={"Accept": "application/json"})
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = result.json()
        self.assertIsNotNone(result_data["genesets"])

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


@parameterized_class(
    [
        {"data_locator": f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad"},
        {"data_locator": f"{FIXTURES_ROOT}/pbmc3k_64.h5ad"},
        {"data_locator": f"{FIXTURES_ROOT}/pbmc3k-CSC-gz.h5ad"},
        {"data_locator": f"{FIXTURES_ROOT}/pbmc3k-CSR-gz.h5ad"},
    ]
)
class EndPointsAnndata(unittest.TestCase, EndPoints):
    """Test Case for endpoints"""

    ANNOTATIONS_ENABLED = False
    GENESETS_READONLY = True

    @classmethod
    def setUpClass(cls):
        if cls == EndPointsAnndata:
            raise unittest.SkipTest("`parameterized_class` bug")

        cls._setupClass(
            cls,
            [
                cls.data_locator,
                "--disable-annotations",
                "--disable-gene-sets-save",
            ],
        )

    @classmethod
    def tearDownClass(cls):
        stop_test_server(cls.ps)

    @property
    def annotations_enabled(self):
        return False

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
        self.assertEqual(len(result_data["positive"]), 7)
        self.assertEqual(len(result_data["negative"]), 7)

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
        self.assertEqual(len(result_data["positive"]), 10)
        self.assertEqual(len(result_data["negative"]), 10)

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


class EndPointsAnndataAnnotations(unittest.TestCase, EndPointsAnnotations):
    """Test Case for endpoints"""

    ANNOTATIONS_ENABLED = True
    GENESETS_READONLY = False

    @classmethod
    def setUpClass(cls):
        cls.data, cls.tmp_dir, cls.annotations = data_with_tmp_annotations(
            MatrixDataType.H5AD, annotations_fixture=True
        )
        cls._setupClass(cls, ["--annotations-file", cls.annotations.label_output_file, cls.data.get_location()])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp_dir)
        stop_test_server(cls.ps)


class EndPointsAnnDataGenesets(unittest.TestCase, EndPoints):
    ANNOTATIONS_ENABLED = False
    GENESETS_READONLY = False

    @classmethod
    def setUpClass(cls):
        cls.tmp_dir = tempfile.mkdtemp()
        genesets_file = path.join(cls.tmp_dir, "test_genesets.csv")
        shutil.copyfile(f"{FIXTURES_ROOT}/pbmc3k-genesets.csv", genesets_file)
        cls._setupClass(
            cls,
            [
                f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad",
                "--disable-annotations",
                "--gene-sets-file",
                genesets_file,
            ],
        )

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp_dir)
        stop_test_server(cls.ps)

    def test_get_genesets_json(self):
        endpoint = "genesets"
        url = f"{self.URL_BASE}{endpoint}"
        result = self.session.get(url, headers={"Accept": "application/json"})
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = result.json()
        self.assertIsNotNone(result_data["genesets"])
        self.assertIsNotNone(result_data["tid"])

        self.assertEqual(
            result_data,
            {
                "genesets": [
                    {
                        "genes": [
                            {"gene_description": " a gene_description", "gene_symbol": "F5"},
                            {"gene_description": "", "gene_symbol": "SUMO3"},
                            {"gene_description": "", "gene_symbol": "SRM"},
                        ],
                        "geneset_description": "a description",
                        "geneset_name": "first gene set name",
                    },
                    {
                        "genes": [
                            {"gene_description": "", "gene_symbol": "RER1"},
                            {"gene_description": "", "gene_symbol": "SIK1"},
                        ],
                        "geneset_description": "",
                        "geneset_name": "second_gene_set",
                    },
                    {"genes": [], "geneset_description": "", "geneset_name": "third gene set"},
                    {"genes": [], "geneset_description": "fourth description", "geneset_name": "fourth_gene_set"},
                    {"genes": [], "geneset_description": "", "geneset_name": "fifth_dataset"},
                    {
                        "genes": [
                            {"gene_description": "", "gene_symbol": "ACD"},
                            {"gene_description": "", "gene_symbol": "AATF"},
                            {"gene_description": "", "gene_symbol": "F5"},
                            {"gene_description": "", "gene_symbol": "PIGU"},
                        ],
                        "geneset_description": "",
                        "geneset_name": "summary test",
                    },
                    {"genes": [], "geneset_description": "", "geneset_name": "geneset_to_delete"},
                    {"genes": [], "geneset_description": "", "geneset_name": "geneset_to_edit"},
                    {
                        "genes": [],
                        "geneset_description": "",
                        "geneset_name": "fill_this_geneset",
                    },
                    {
                        "genes": [{"gene_description": "", "gene_symbol": "SIK1"}],
                        "geneset_description": "",
                        "geneset_name": "empty_this_geneset",
                    },
                    {
                        "genes": [{"gene_description": "", "gene_symbol": "SIK1"}],
                        "geneset_description": "",
                        "geneset_name": "brush_this_gene",
                    },
                ],
                "tid": 0,
            },
        )

    def test_get_genesets_csv(self):
        endpoint = "genesets"
        url = f"{self.URL_BASE}{endpoint}"
        result = self.session.get(url, headers={"Accept": "text/csv"})
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "text/csv")
        self.assertEqual(
            result.text,
            """gene_set_name,gene_set_description,gene_symbol,gene_description\r
first gene set name,a description,F5, a gene_description\r
first gene set name,a description,SUMO3,\r
first gene set name,a description,SRM,\r
second_gene_set,,RER1,\r
second_gene_set,,SIK1,\r
third gene set,,,\r
fourth_gene_set,fourth description,,\r
fifth_dataset,,,\r
summary test,,ACD,\r
summary test,,AATF,\r
summary test,,F5,\r
summary test,,PIGU,\r
geneset_to_delete,,,\r
geneset_to_edit,,,\r
fill_this_geneset,,,\r
empty_this_geneset,,SIK1,\r
brush_this_gene,,SIK1,\r
""",
        )

    def test_put_genesets(self):
        endpoint = "genesets"
        url = f"{self.URL_BASE}{endpoint}"

        # assume we start with TID 0
        result = self.session.get(url, headers={"Accept": "application/json"})
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.json()["tid"], 0)

        test1 = {"tid": 3, "genesets": []}
        result = self.session.put(url, json=test1)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        result = self.session.get(url, headers={"Accept": "application/json"})
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.json(), test1)

        # stale TID
        result = self.session.put(url, json=test1)
        self.assertEqual(result.status_code, HTTPStatus.NOT_FOUND)

        test2 = {
            "tid": 4,
            "genesets": [
                {"geneset_name": "foobar", "genes": []},
                {"geneset_name": "contains a space", "genes": []},
                {"geneset_name": "contains_weird_characters: #$%^&*()_+=-!@<>,./?';:\"[]{}|\\", "genes": []},
            ],
        }
        test2_response = {
            "tid": 4,
            "genesets": [
                {"geneset_name": "foobar", "geneset_description": "", "genes": []},
                {"geneset_name": "contains a space", "geneset_description": "", "genes": []},
                {
                    "geneset_name": "contains_weird_characters: #$%^&*()_+=-!@<>,./?';:\"[]{}|\\",
                    "geneset_description": "",
                    "genes": [],
                },
            ],
        }
        result = self.session.put(url, json=test2)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        result = self.session.get(url, headers={"Accept": "application/json"})
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.json(), test2_response)

        test3 = {
            "tid": 5,
            "genesets": [
                {
                    "geneset_name": "foobar",
                    "geneset_description": "",
                    "genes": [
                        {
                            "gene_symbol": "F5",
                            "gene_description": "",
                        }
                    ],
                }
            ],
        }
        result = self.session.put(url, json=test3)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        result = self.session.get(url, headers={"Accept": "application/json"})
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.json(), test3)

    def test_put_genesets_malformed(self):
        """test malformed submissions that we expect the backend to catch/tolerate"""
        endpoint = "genesets"
        url = f"{self.URL_BASE}{endpoint}"

        result = self.session.get(url, headers={"Accept": "application/json"})
        self.assertEqual(result.status_code, HTTPStatus.OK)
        original_data = result.json()
        tid = original_data["tid"]

        def test_case(test, expected_code, original_data):
            """check for expected error AND that no change was made to the original state"""
            result = self.session.put(url, json=test)
            self.assertEqual(result.status_code, expected_code)
            result = self.session.get(url, headers={"Accept": "application/json"})
            self.assertEqual(result.status_code, HTTPStatus.OK)
            self.assertEqual(result.json(), original_data)

        # missing or malformed genesets
        test_case(
            {"tid": tid + 1},
            HTTPStatus.BAD_REQUEST,
            original_data,
        )
        test_case(
            {"tid": tid + 1, "genesets": 99},
            HTTPStatus.BAD_REQUEST,
            original_data,
        )

        # illegal geneset_name
        test_case(
            {"tid": tid + 1, "genesets": [{"geneset_name": " foo", "genes": []}]},
            HTTPStatus.BAD_REQUEST,
            original_data,
        )
        test_case(
            {"tid": tid + 1, "genesets": [{"geneset_name": "foo ", "genes": []}]},
            HTTPStatus.BAD_REQUEST,
            original_data,
        )
        test_case(
            {"tid": tid + 1, "genesets": [{"geneset_name": "f  oo", "genes": []}]},
            HTTPStatus.BAD_REQUEST,
            original_data,
        )
        test_case(
            {"tid": tid + 1, "genesets": [{"geneset_name": "f\too", "genes": []}]},
            HTTPStatus.BAD_REQUEST,
            original_data,
        )
        test_case(
            {"tid": tid + 1, "genesets": [{"geneset_name": "f\roo", "genes": []}]},
            HTTPStatus.BAD_REQUEST,
            original_data,
        )
        test_case(
            {"tid": tid + 1, "genesets": [{"geneset_name": "f\noo", "genes": []}]},
            HTTPStatus.BAD_REQUEST,
            original_data,
        )
        test_case(
            {"tid": tid + 1, "genesets": [{"geneset_name": "f\voo", "genes": []}]},
            HTTPStatus.BAD_REQUEST,
            original_data,
        )

        # duplicate geneset_name
        test_case(
            {
                "tid": tid + 1,
                "genesets": [
                    {"geneset_name": "foo", "genes": []},
                    {"geneset_name": "foo", "genes": []},
                ],
            },
            HTTPStatus.BAD_REQUEST,
            original_data,
        )

        # missing geneset_name
        test_case(
            {"tid": tid + 1, "genesets": [{"genes": []}]},
            HTTPStatus.BAD_REQUEST,
            original_data,
        )

        # non-numeric TID
        test_case(
            {"tid": [], "genesets": [{"geneset_name": "foo", "genes": []}]},
            HTTPStatus.BAD_REQUEST,
            original_data,
        )
        test_case(
            {"tid": None, "genesets": [{"geneset_name": "foo", "genes": []}]},
            HTTPStatus.BAD_REQUEST,
            original_data,
        )
        test_case(
            {"tid": "not a number", "genesets": [{"geneset_name": "foo", "genes": []}]},
            HTTPStatus.BAD_REQUEST,
            original_data,
        )

        # duplicate gene_symbol
        test_case(
            {
                "tid": "not a number",
                "genesets": [{"geneset_name": "foo", "genes": [{"gene_symbol": "SIK1"}, {"gene_symbol": "SIK1"}]}],
            },
            HTTPStatus.BAD_REQUEST,
            original_data,
        )

        # gene_symbol is not a string
        test_case(
            {
                "tid": "not a number",
                "genesets": [{"geneset_name": "foo", "genes": [{"gene_symbol": 99}]}],
            },
            HTTPStatus.BAD_REQUEST,
            original_data,
        )

    def test_get_geneset_summary_edge_cases(self):
        # attempt to summarize _all_ genesets, including edge cases with zero or one gene
        result = self.session.get(f"{self.URL_BASE}genesets", headers={"Accept": "application/json"})
        self.assertEqual(result.status_code, HTTPStatus.OK)
        genesets = result.json()["genesets"]

        endpoint = "summarize/var"
        index_col_name = self.schema["schema"]["annotations"]["var"]["index"]
        for gs in genesets:
            genes = [g["gene_symbol"] for g in gs["genes"]]
            filter = "&".join([f"var:{index_col_name}={gene}" for gene in genes])
            query = f"method=mean&{filter}"
            query_hash = hashlib.sha1(query.encode()).hexdigest()
            url = f"{self.URL_BASE}{endpoint}?{query}"

            result = self.session.get(url, headers={"Accept": "application/octet-stream"})
            self.assertEqual(result.status_code, HTTPStatus.OK)
            self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
            df = decode_fbs.decode_matrix_FBS(result.content)
            self.assertEqual(df["n_rows"], 2638)
            self.assertEqual(df["n_cols"], 1)
            self.assertEqual(df["col_idx"], [query_hash])
