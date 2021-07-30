import json
import os
import time
from http import HTTPStatus
import hashlib

import requests


from backend.czi_hosted.common.config.app_config import AppConfig
from backend.test import decode_fbs
from backend.test.fixtures.fixtures import pbmc3k_colors
from backend.test.test_czi_hosted.unit import BaseTest, skip_if

BAD_FILTER = {"filter": {"obs": {"annotation_value": [{"name": "xyz"}]}}}

class EndPoints(BaseTest):
    @classmethod
    def setUpClass(cls, app_config=None):
        super().setUpClass(app_config)
        cls.app.testing = True
        cls.client = cls.app.test_client()
        os.environ["SKIP_STATIC"] = "True"
        for i in range(90):
            try:
                result = cls.client.get(f"{cls.TEST_URL_BASE}schema")
                cls.schema = json.loads(result.data)
            except requests.exceptions.ConnectionError:
                time.sleep(1)

    def test_initialize(self):
        endpoint = "schema"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        result = self.client.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = json.loads(result.data)
        self.assertEqual(result_data["schema"]["dataframe"]["nObs"], 2638)
        self.assertEqual(len(result_data["schema"]["annotations"]["obs"]), 2)
        self.assertEqual(
            len(result_data["schema"]["annotations"]["obs"]["columns"]), 5
        )

    def test_config(self):
        endpoint = "config"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        result = self.client.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = json.loads(result.data)
        self.assertIn("library_versions", result_data["config"])
        self.assertEqual(result_data["config"]["displayNames"]["dataset"], "pbmc3k")

    def test_get_layout_fbs(self):
        endpoint = "layout/obs"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 8)
        self.assertIsNotNone(df["columns"])
        self.assertSetEqual(
            set(df["col_idx"]),
            {"pca_0", "pca_1", "tsne_0", "tsne_1", "umap_0", "umap_1", "draw_graph_fr_0", "draw_graph_fr_1"},
        )
        self.assertIsNone(df["row_idx"])
        self.assertEqual(len(df["columns"]), df["n_cols"])

    def test_bad_filter(self):
        endpoint = "data/var"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.put(url, headers=header, json=BAD_FILTER)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

    def test_get_annotations_obs_fbs(self):
        endpoint = "annotations/obs"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 5)
        self.assertIsNotNone(df["columns"])
        self.assertIsNone(df["row_idx"])
        self.assertEqual(len(df["columns"]), df["n_cols"])
        obs_index_col_name = self.schema["schema"]["annotations"]["obs"]["index"]
        self.assertCountEqual(
            df["col_idx"],
            [obs_index_col_name, "n_genes", "percent_mito", "n_counts", "louvain"],
        )

    def test_get_annotations_obs_keys_fbs(self):
        endpoint = "annotations/obs"
        query = "annotation-name=n_genes&annotation-name=percent_mito"
        url = f"{self.TEST_URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 2)
        self.assertIsNotNone(df["columns"])
        self.assertIsNone(df["row_idx"])
        self.assertEqual(len(df["columns"]), df["n_cols"])
        self.assertCountEqual(df["col_idx"], ["n_genes", "percent_mito"])

    def test_get_annotations_obs_error(self):
        endpoint = "annotations/obs"
        query = "annotation-name=notakey"
        url = f"{self.TEST_URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

# TEMP: Testing count 15 to match hardcoded values for diffexp
# TODO(#1281): Switch back to dynamic values
    def test_diff_exp(self):
        endpoint = "diffexp/obs"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        params = {
            "mode": "topN",
            "set1": {"filter": {"obs": {"annotation_value": [{"name": "louvain", "values": ["NK cells"]}]}}},
            "set2": {"filter": {"obs": {"annotation_value": [{"name": "louvain", "values": ["CD8 T cells"]}]}}},
            "count": 15,
        }
        result = self.client.post(url, json=params)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = json.loads(result.data)
        self.assertEqual(len(result_data['positive']), 15)
        self.assertEqual(len(result_data['negative']), 15)

    def test_diff_exp_indices(self):
        endpoint = "diffexp/obs"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        params = {
            "mode": "topN",
            "count": 15,
            "set1": {"filter": {"obs": {"index": [[0, 500]]}}},
            "set2": {"filter": {"obs": {"index": [[500, 1000]]}}},
        }
        result = self.client.post(url, json=params)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = json.loads(result.data)
        self.assertEqual(len(result_data['positive']), 15)
        self.assertEqual(len(result_data['negative']), 15)

    def test_get_annotations_var_fbs(self):
        endpoint = "annotations/var"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
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
        url = f"{self.TEST_URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertEqual(df["n_rows"], 1838)
        self.assertEqual(df["n_cols"], 1)
        self.assertIsNotNone(df["columns"])
        self.assertIsNone(df["row_idx"])
        self.assertEqual(len(df["columns"]), df["n_cols"])
        self.assertCountEqual(df["col_idx"], ["n_cells"])

    def test_get_annotations_var_error(self):
        endpoint = "annotations/var"
        query = "annotation-name=notakey"
        url = f"{self.TEST_URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

    def test_data_mimetype_error(self):
        endpoint = "data/var"
        header = {"Accept": "xxx"}
        url = f"{self.TEST_URL_BASE}{endpoint}"
        result = self.client.put(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.NOT_ACCEPTABLE)

    def test_fbs_default(self):
        endpoint = "data/var"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        headers = {"Accept": "application/octet-stream"}
        result = self.client.put(url, headers=headers)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

        filter = {"filter": {"var": {"index": [0, 1, 4]}}}
        result = self.client.put(url, headers=headers, json=filter)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")

    def test_data_put_fbs(self):
        endpoint = "data/var"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.put(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

    def test_data_get_fbs(self):
        endpoint = "data/var"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.BAD_REQUEST)

    def test_data_put_filter_fbs(self):
        endpoint = "data/var"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        filter = {"filter": {"var": {"index": [0, 1, 4]}}}
        result = self.client.put(url, headers=header, json=filter)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
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
        url = f"{self.TEST_URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 1)

    def test_data_get_unknown_filter_fbs(self):
        index_col_name = self.schema["schema"]["annotations"]["var"]["index"]
        endpoint = "data/var"
        query = f"var:{index_col_name}=UNKNOWN"
        url = f"{self.TEST_URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 0)

    def test_data_put_single_var(self):
        endpoint = "data/var"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        index_col_name = self.schema["schema"]["annotations"]["var"]["index"]
        var_filter = {"filter": {"var": {"annotation_value": [{"name": index_col_name, "values": ["RER1"]}]}}}
        result = self.client.put(url, headers=header, json=var_filter)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 1)

    def test_colors(self):
        endpoint = "colors"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        result = self.client.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = json.loads(result.data)
        self.assertEqual(result_data, pbmc3k_colors)

    @skip_if(lambda x: os.getenv("SKIP_STATIC"), "Skip static test when running locally")
    def test_static(self):
        endpoint = "static"
        file = "assets/favicon.ico"
        url = f"{endpoint}/{file}"
        result = self.client.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)

    def test_genesets_config(self):
        result = self.client.get(f"{self.TEST_URL_BASE}config")
        config_data = json.loads(result.data)
        params = config_data["config"]["parameters"]
        annotations_genesets = params["annotations_genesets"]
        annotations_genesets_readonly = params["annotations_genesets_readonly"]
        annotations_genesets_summary_methods = params["annotations_genesets_summary_methods"]
        self.assertTrue(annotations_genesets)
        self.assertTrue(annotations_genesets_readonly)
        self.assertEqual(annotations_genesets_summary_methods, ["mean"])

    def test_get_genesets(self):
        endpoint = "genesets"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        result = self.client.get(url, headers={"Accept": "application/json"})
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = json.loads(result.data)
        self.assertIsNotNone(result_data["genesets"])

    def test_get_summaryvar(self):
        index_col_name = self.schema["schema"]["annotations"]["var"]["index"]
        endpoint = "summarize/var"

        # single column
        filter = f"var:{index_col_name}=F5"
        query = f"method=mean&{filter}"
        query_hash = hashlib.sha1(query.encode()).hexdigest()
        url = f"{self.TEST_URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 1)
        self.assertEqual(df["col_idx"], [query_hash])
        self.assertAlmostEqual(df["columns"][0][0], -0.110451095)

        # multi-column
        col_names = ["F5", "BEB3", "SIK1"]
        filter = "&".join([f"var:{index_col_name}={name}" for name in col_names])
        query = f"method=mean&{filter}"
        query_hash = hashlib.sha1(query.encode()).hexdigest()
        url = f"{self.TEST_URL_BASE}{endpoint}?{query}"
        header = {"Accept": "application/octet-stream"}
        result = self.client.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
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
        url = f"{self.TEST_URL_BASE}{endpoint}?key={query_hash}"
        result = self.client.post(url, headers=headers, data=query)

        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 1)
        self.assertEqual(df["col_idx"], [query_hash])
        self.assertAlmostEqual(df["columns"][0][0], -0.110451095)

        # multi-column
        col_names = ["F5", "BEB3", "SIK1"]
        filter = "&".join([f"var:{index_col_name}={name}" for name in col_names])
        query = f"method=mean&{filter}"
        query_hash = hashlib.sha1(query.encode()).hexdigest()
        url = f"{self.TEST_URL_BASE}{endpoint}?key={query_hash}"
        result = self.client.post(url, headers=headers, data=query)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertEqual(df["n_rows"], 2638)
        self.assertEqual(df["n_cols"], 1)
        self.assertEqual(df["col_idx"], [query_hash])
        self.assertAlmostEqual(df["columns"][0][0], -0.16628358)



class EndPointsCxg(EndPoints):
    """Test Case for endpoints"""
    @classmethod
    def setUpClass(cls):
        app_config = AppConfig()
        app_config.update_default_dataset_config(user_annotations__enable=False)

    def test_get_genesets_json(self):
        self.app.auth.is_user_authenticated = lambda: True
        endpoint = "genesets"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        result = self.client.get(url, headers={"Accept": "application/json"})
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")
        result_data = json.loads(result.data)
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
                    {'genes': [], 'geneset_description': '', 'geneset_name': 'geneset_to_delete'},
                    {'genes': [], 'geneset_description': '', 'geneset_name': 'geneset_to_edit'},
                    {
                        'genes': [],
                        'geneset_description': '',
                        'geneset_name': 'fill_this_geneset'
                    },
                    {
                        'genes': [{'gene_description': '', 'gene_symbol': 'SIK1'}],
                        'geneset_description': '',
                        'geneset_name': 'empty_this_geneset'
                    },
                    {
                        'genes': [{'gene_description': '', 'gene_symbol': 'SIK1'}],
                        'geneset_description': '',
                        'geneset_name': 'brush_this_gene'
                    }
                ],
                "tid": 0,
            },
        )

    def test_get_genesets_csv(self):
        endpoint = "genesets"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        self.app.auth.is_user_authenticated = lambda: True
        result = self.client.get(url, headers={"Accept": "text/csv"})
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "text/csv")
        expected_data = """gene_set_name,gene_set_description,gene_symbol,gene_description\r
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
"""
        self.assertEqual(result.data.decode("utf-8"), expected_data)

    def test_put_genesets(self):
        endpoint = "genesets"
        url = f"{self.TEST_URL_BASE}{endpoint}"

        result = self.client.get(url, headers={"Accept": "application/json"})
        self.assertEqual(result.status_code, HTTPStatus.OK)

        test1 = {"tid": 3, "genesets": []}
        result = self.client.put(url, json=test1)

        self.assertEqual(result.status_code, HTTPStatus.METHOD_NOT_ALLOWED)

