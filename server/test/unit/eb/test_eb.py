import unittest
import tempfile
import requests
import subprocess
from server.test import PROJECT_ROOT, FIXTURES_ROOT
from server.common.config.app_config import AppConfig
from contextlib import contextmanager
import time


@contextmanager
def run_eb_app(tempdirname):
    ps = subprocess.Popen(["python", "artifact.dir/application.py"], cwd=tempdirname)
    server = "http://localhost:5000"
    for _ in range(10):
        try:
            requests.get(f"{server}/health")
            break
        except requests.exceptions.ConnectionError:
            time.sleep(1)

    try:
        yield server
    finally:
        try:
            ps.terminate()
        except ProcessLookupError:
            pass


class Elastic_Beanstalk_Test(unittest.TestCase):
    def test_run(self):

        tempdir = tempfile.TemporaryDirectory(dir=f"{PROJECT_ROOT}/server")
        tempdirname = tempdir.name

        c = AppConfig()
        # test that eb works
        c.update_server_config(
            multi_dataset__dataroot=f"{FIXTURES_ROOT}", app__flask_secret_key="open sesame"
        )

        c.complete_config()
        c.write_config(f"{tempdirname}/config.yaml")

        subprocess.check_call(f"git ls-files . | cpio -pdm {tempdirname}", cwd=f"{PROJECT_ROOT}/server/eb", shell=True)
        subprocess.check_call(["make", "build"], cwd=tempdirname)

        with run_eb_app(tempdirname) as server:
            session = requests.Session()

            r = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/config")
            data_config = r.json()
            assert data_config["config"]["displayNames"]["dataset"] == "pbmc3k"
