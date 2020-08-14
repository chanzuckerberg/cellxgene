import unittest
from multiprocessing import Process, Queue

from click.testing import CliRunner

from server.cli.launch import launch


class TestLaunch(unittest.TestCase):

    def run_cli_command(self, cli_args, q):
        runner = CliRunner()
        q.put(runner.invoke(launch, cli_args))

    def test__launch__simple(self):
        q = Queue()
        p = Process(target=self.run_cli_command, name="launch_cellxgene",
                    args=(['example-dataset/pbmc3k.h5ad'], q,))
        p.start()
        import time
        time.sleep(10)
        runner_result = q.get()
        p.terminate()
        self.assertIn('Launching!', runner_result.output)
        assert runner_result.exit_code == 0
