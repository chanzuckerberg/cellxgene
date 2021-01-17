import unittest
import tempfile
import requests
import subprocess
from local_server.test import PROJECT_ROOT, FIXTURES_ROOT
from local_server.common.config.app_config import AppConfig
from contextlib import contextmanager
import time
import os


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

        tempdir = tempfile.TemporaryDirectory(dir=f"{PROJECT_ROOT}/local_server")
        tempdirname = tempdir.name

        config = AppConfig()
        # test that eb works
        config.update_server_config(multi_dataset__dataroot=f"{FIXTURES_ROOT}", app__flask_secret_key="open sesame")

        config.complete_config()
        config.write_config(f"{tempdirname}/config.yaml")

        subprocess.check_call(f"git ls-files . | cpio -pdm {tempdirname}", cwd=f"{PROJECT_ROOT}/local_server/eb", shell=True)
        subprocess.check_call(["make", "build"], cwd=tempdirname)

        with run_eb_app(tempdirname) as server:
            session = requests.Session()

            r = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/config")
            data_config = r.json()
            assert data_config["config"]["displayNames"]["dataset"] == "pbmc3k"

    def test_config(self):
        check_config_script = os.path.join(PROJECT_ROOT, "local_server", "eb", "check_config.py")
        with tempfile.TemporaryDirectory() as tempdir:
            configfile = os.path.join(tempdir, "config.yaml")
            app_config = AppConfig()
            app_config.update_server_config(multi_dataset__dataroot=f"{FIXTURES_ROOT}")
            app_config.write_config(configfile)

            command = ["python", check_config_script, configfile]

            # test failure mode (flask_secret_key not set)
            env = os.environ.copy()
            env.pop("CXG_SECRET_KEY", None)
            with self.assertRaises(subprocess.CalledProcessError) as exception_context:
                subprocess.check_output(command, env=env)
            output = str(exception_context.exception.stdout, "utf-8")
            self.assertTrue(
                output.startswith(
                    "Error: Invalid type for attribute: app__flask_secret_key, expected type str, got NoneType"
                )
            )
            self.assertEqual(exception_context.exception.returncode, 1)

            # test passing case
            env = os.environ.copy()
            env["CXG_SECRET_KEY"] = "secret"
            output = subprocess.check_output(command, env=env)
            output = str(output, "utf-8")
            self.assertTrue(output.startswith("PASS"))
