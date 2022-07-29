import os
import shutil
import unittest
from tempfile import mkstemp, TemporaryDirectory

import mlflow
from click.testing import CliRunner

from server.cli.annotate import annotate
from test.unit.cli.fixtures.mlflow_model_fixture import FakeModel


def write_model(model) -> str:
    with TemporaryDirectory() as mlflow_model_dir:
        fixtures_path = os.path.join(os.path.dirname(__file__), 'fixtures')
        mlflow.pyfunc.save_model(mlflow_model_dir,
                                 loader_module='fixtures',
                                 code_path=[fixtures_path])
        return shutil.make_archive(mkstemp()[1], "zip", mlflow_model_dir)


class TestCliAnnotate(unittest.TestCase):
    def test__annotate__loads_and_runs(self):
        """
        Invokes the `annotate` subcommand of cellxgene CLI, using a CliRunner() programmatic invocation.

        This tests the happy path case:
        1) Command line options are parsed;
        2) An MLflow model zip archive can be read in (from local disk), unpacked, and invoked;
        3) The correct options are passed to the MLflow model.
        4) The annotate subcommand exits successfully.

        This does not verify model output or predictions (it's a fake MLflow model, after all); it's up to the real model
        to output its predictions as it wants, but this is specific to the model and so not tested here.

        The  CliRunner() invokes the subcommand in a subprocess, and the annotate subcommand itself invokes the MLflow
        model in yet another subprocess. So while this test can help determine if everything is working, it is not a
        simple matter to debug in the case of a failure. However, the stdout/stderr of the MLflow process is captured
        by the CliRunner() subprocess, so errors can be inspected in result.stdout when debugging this test. Hope this
        helps!
        """

        _, query_dataset_file_path = mkstemp()
        model_file_path = write_model(FakeModel())

        result = CliRunner().invoke(
            annotate,
            [
                "--input-h5ad-file",
                query_dataset_file_path,
                "--model-url",
                model_file_path,
                "--output-h5ad-file",
                f"{query_dataset_file_path}.output",
                # avoid having mflow create conda env or virtualenv when in test env;
                # this avoids making pip remote requests and is also faster
                "--mlflow-env-manager", "local"
            ],
        )

        # to help debugging, show the output from the CliRunner and MLflow stdout
        if result.exit_code:
            print(result.stdout)

        self.assertEqual(0, result.exit_code, "runs successfully")

        # The FakeModel will print it inputs to stdout, as "__MODEL_INPUT__={...}", allowing us to assert that it received valid inputs.
        self.assertIn(
            "__MODEL_INPUT__={"
            f'"query_dataset_h5ad_path": "{query_dataset_file_path}", '
            f'"output_h5ad_path": "{query_dataset_file_path}.output", '
            '"annotation_prefix": "cxg_cell_type", "classifier": "default", '
            '"organism": "Homo sapiens", "use_gpu": true}',
            result.stdout,
            "inputs passed correctly",
        )

    def test__annotate__verifies_mutually_exclusive_options(self):
        required_options = ["--input-h5ad-file", "some.h5ad", "--model-url", "some_url"]
        result = CliRunner().invoke(
            annotate,
            required_options + [],
        )

        self.assertNotEqual(0, result.exit_code, "aborts with non-success code")
        self.assertIn(
            "--update_h5ad_file or --output_h5ad_file must be specified",
            result.stdout,
            "error message displayed",
        )

        result = CliRunner().invoke(
            annotate, required_options + ["--output-h5ad-file", "some_arg", "--update-h5ad-file"]
        )

        self.assertNotEqual(0, result.exit_code, "aborts with non-success code")
        self.assertIn(
            "--update_h5ad_file and --output_h5ad_file are mutually exclusive",
            result.stdout,
            "error message displayed",
        )


# TODO:
# Test annotate cli args more comprehensively
# Test server.cli.annotate._validate_options
# Test model caching feature works
# Test model loading from s3 works (maybe w/just a real model)

if __name__ == "__main__":
    unittest.main()
