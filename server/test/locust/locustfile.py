import json
import random

import requests
from config import DataSets
from locust import HttpUser, SequentialTaskSet, task, between, TaskSet
from locust.clients import HttpSession
from requests.packages.urllib3.exceptions import InsecureRequestWarning

import server.test.unit.decode_fbs as decode_fbs

requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

"""
Simple locust stress test definition for cellxgene
"""

API_SUFFIX = "api/v0.2"


class CellXGeneTasks(TaskSet):
    """
    Simulate use against a single dataset
    """

    def on_start(self):

        self.client.verify = False
        self.dataset = random.choice(DataSets)

        with self.client.get(
            f"{self.dataset}/{API_SUFFIX}/schema", stream=True, catch_response=True
        ) as schema_response:
            if schema_response.status_code == 200:
                self.schema = schema_response.json()["schema"]
            else:
                self.schema = None

        with self.client.get(
            f"{self.dataset}/{API_SUFFIX}/config", stream=True, catch_response=True
        ) as config_response:
            if config_response.status_code == 200:
                self.config = config_response.json()["config"]
            else:
                self.config = None

        with self.client.get(
            f"{self.dataset}/{API_SUFFIX}/annotations/var?annotation-name={self.var_index_name()}",
            headers={"Accept": "application/octet-stream"},
            catch_response=True,
        ) as var_index_response:
            if var_index_response.status_code == 200:
                df = decode_fbs.decode_matrix_FBS(var_index_response.content)
                gene_names_idx = df["col_idx"].index(self.var_index_name())
                self.gene_names = df["columns"][gene_names_idx]
            else:
                self.gene_names = []

    def var_index_name(self):
        if self.schema is not None:
            return self.schema["annotations"]["var"]["index"]
        return None

    def obs_annotation_names(self):
        if self.schema is not None:
            return [col["name"] for col in self.schema["annotations"]["obs"]["columns"]]
        return []

    def layout_names(self):
        if self.schema is not None:
            return [layout["name"] for layout in self.schema["layout"]["obs"]]
        else:
            return []

    @task(2)
    class InitializeClient(SequentialTaskSet):
        """
        Initial loading of cellxgene - when the user hits the main route.

        Currently this sequence skips some of the static assets, which are quite small and should be served by the
        HTTP server directly.

        1. Load index.html, etc.
        2. Concurrently load /config, /schema
        3. Concurrently load /layout/obs, /annotations/var?annotation-name=<the index>
        -- Does initial render --
        4. Concurrently load all /annotations/obs and all /layouts/obs
        -- Fully initialized --
        """

        # Users hit all of the init routes as fast as they can, subject to the ordering constraints and network latency.
        wait_time = between(0.01, 0.1)

        def on_start(self):
            self.dataset = self.parent.dataset
            self.client.verify = False
            self.api_less_client = HttpSession(
                base_url=self.client.base_url.replace("api.", "").replace("cellxgene/", ""),
                request_success=self.client.request_success,
                request_failure=self.client.request_failure,
            )

        @task
        def index(self):
            self.api_less_client.get(f"{self.dataset}", stream=True)

        @task
        def loadConfigAndSchema(self):
            self.client.get(f"{self.dataset}/{API_SUFFIX}/schema", stream=True, catch_response=True)
            self.client.get(f"{self.dataset}/{API_SUFFIX}/config", stream=True, catch_response=True)

        @task
        def loadBootstrapData(self):
            self.client.get(
                f"{self.dataset}/{API_SUFFIX}/layout/obs", headers={"Accept": "application/octet-stream"}, stream=True
            )
            self.client.get(
                f"{self.dataset}/{API_SUFFIX}/annotations/var?annotation-name={self.parent.var_index_name()}",
                headers={"Accept": "application/octet-stream"},
                catch_response=True,
            )

        @task
        def loadObsAnnotationsAndLayouts(self):
            obs_names = self.parent.obs_annotation_names()
            for name in obs_names:
                self.client.get(
                    f"{self.dataset}/{API_SUFFIX}/annotations/obs?annotation-name={name}",
                    headers={"Accept": "application/octet-stream"},
                    stream=True,
                )

            layouts = self.parent.layolut_names()
            for name in layouts:
                self.client.get(
                    f"{self.dataset}/{API_SUFFIX}/annotations/obs?layout-name={name}",
                    headers={"Accept": "application/octet-stream"},
                    stream=True,
                )

        @task
        def done(self):
            self.interrupt()

    @task(1)
    def load_expression(self):
        """
        Simulate user occasionally loading some expression data for a gene
        """

        gene_name = random.choice(self.gene_names)
        filter = {"filter": {"var": {"annotation_value": [{"name": self.var_index_name(), "values": [gene_name]}]}}}
        self.client.put(
            f"{self.dataset}/{API_SUFFIX}/data/var",
            data=json.dumps(filter),
            headers={"Content-Type": "application/json", "Accept": "application/octet-stream"},
            stream=True,
        ).close()


class CellxgeneUser(HttpUser):
    tasks = [CellXGeneTasks]

    # Most ops do not require back-end interaction, so slow cadence for users
    wait_time = between(10, 60)
