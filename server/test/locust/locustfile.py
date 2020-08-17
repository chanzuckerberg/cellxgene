from locust import HttpLocust, TaskSet, TaskSequence, seq_task, task
from locust.wait_time import between
import random
import json
from gevent.pool import Group

import server.test.unit.decode_fbs as decode_fbs
from config import DataSets


"""
Simple locust stress test defition for cellxgene
"""
API = "/api/v0.2"


class ViewDataset(TaskSet):
    """
    Simulate use against a single dataset
    """

    def on_start(self):
        self.client.verify = False
        self.dataset = random.choice(DataSets)

        with self.client.get(f"{self.dataset}{API}/config", catch_response=True) as r:
            if r.status_code == 200:
                self.config = r.json()["config"]
                r.success()
            else:
                self.config = None
                r.failure(f"bad response code {r.status_code}")

        with self.client.get(f"{self.dataset}{API}/schema", catch_response=True) as r:
            if r.status_code == 200:
                self.schema = r.json()["schema"]
                r.success()
            else:
                self.schema = None
                r.failure(f"bad response code {r.status_code}")

        with self.client.get(
            f"{self.dataset}{API}/annotations/var?annotation-name={self.var_index_name()}",
            headers={"Accept": "application/octet-stream"},
            catch_response=True,
        ) as r:
            if r.status_code == 200:
                df = decode_fbs.decode_matrix_FBS(r.content)
                gene_names_idx = df["col_idx"].index(self.var_index_name())
                self.gene_names = df["columns"][gene_names_idx]
            else:
                self.gene_names = None
                r.failure(f"bad response code {r.status_code}")

    def var_index_name(self):
        if self.schema is None:
            return None
        return self.schema["annotations"]["var"]["index"]

    def obs_annotation_names(self):
        if self.schema is None:
            return []
        return [col["name"] for col in self.schema["annotations"]["obs"]["columns"]]

    @task(2)
    class InitializeClient(TaskSequence):
        """
        Initial loading of cellxgene - when the user hits the main route.

        Currently this sequence skips some of the static assets, which are quite
        small and should be served by the HTTP server directly.

        1. load index.html, etc.
        2. concurrently load /config, /schema
        3. concurrently load /layout/obs, /annotations/var?annotation-name=<the index>
        -- does intitial render --
        4. concurrently load all /annotations/obs
        -- fully initialized --
        """

        # users hit all of the init routes as fast as they can, subject to the ordering constraints
        # and network latency
        wait_time = between(0.01, 0.1)

        def on_start(self):
            self.dataset = self.parent.dataset
            self.client.verify = False

        @seq_task(1)
        def index(self):
            self.client.get(f"{self.dataset}/", stream=True).close()

        @seq_task(2)
        def loadConfigSchema(self):
            def config():
                self.client.get(f"{self.dataset}{API}/config", stream=True).close()

            def schema():
                self.client.get(f"{self.dataset}{API}/schema", stream=True).close()

            group = Group()
            group.spawn(config)
            group.spawn(schema)
            group.join()

        @seq_task(3)
        def loadBootstrapData(self):
            def layout():
                self.client.get(
                    f"{self.dataset}{API}/layout/obs", headers={"Accept": "application/octet-stream"}, stream=True
                ).close()

            def varAnnotationIndex():
                self.client.get(
                    f"{self.dataset}{API}/annotations/var?annotation-name={self.parent.var_index_name()}",
                    headers={"Accept": "application/octet-stream"},
                    stream=True,
                ).close()

            group = Group()
            group.spawn(layout)
            group.spawn(varAnnotationIndex)
            group.join()

        @seq_task(4)
        def loadObsAnnotations(self):
            def obs_annotation(name):
                self.client.get(
                    f"{self.dataset}{API}/annotations/obs?annotation-name={name}",
                    headers={"Accept": "application/octet-stream"},
                    stream=True,
                ).close()

            obs_names = self.parent.obs_annotation_names()
            group = Group()
            for name in obs_names:
                group.spawn(obs_annotation, name)
            group.join()

        @seq_task(5)
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
            f"{self.dataset}{API}/data/var",
            data=json.dumps(filter),
            headers={"Content-Type": "application/json", "Accept": "application/octet-stream"},
            stream=True,
        ).close()


class CellxgeneUser(HttpLocust):
    task_set = ViewDataset

    # most ops do not require back-end interaction, so slow cadence
    # for users
    wait_time = between(10, 60)
