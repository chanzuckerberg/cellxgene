import os
import shutil
import unittest
from tempfile import mkstemp, TemporaryDirectory, NamedTemporaryFile

import mlflow
from click.testing import CliRunner

from server.cli.annotate import annotate
from test.unit.cli.fixtures.mlflow_model_fixture import FakeModel


def write_model(model) -> str:
    with TemporaryDirectory() as mlflow_model_dir:
        fixtures_path = os.path.join(os.path.dirname(__file__), "fixtures")
        mlflow.pyfunc.save_model(mlflow_model_dir, loader_module="fixtures", code_path=[fixtures_path])
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
                query_dataset_file_path,
                "--model-url",
                model_file_path,
                "--output-h5ad-file",
                f"{query_dataset_file_path}.output",
                # avoid having mflow create conda env or virtualenv when in test env;
                # this avoids making pip remote requests and is also faster
                "--mlflow-env-manager",
                "local",
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
        self.assertIn(
            f"Wrote annotations to {query_dataset_file_path}.output",
            result.stdout,
            "success message is correct",
        )

    def test__annotate__requires_overwrite_option_when_output_file_exists(self):

        with NamedTemporaryFile() as input_h5ad, NamedTemporaryFile() as existing_file:
            required_options = [input_h5ad.name, "--output-h5ad-file", existing_file.name, "--model-url", "some_url"]
            result = CliRunner().invoke(
                annotate,
                required_options + [],
            )

            self.assertNotEqual(0, result.exit_code, "aborts with non-success code")
            self.assertIn(
                "try using the flag --overwrite",
                result.stdout,
                "error message displayed",
            )

    def test__annotate__overwrite_option_allows_overwrite_of_existing_output_file(self):
        model_file_path = write_model(FakeModel())

        with NamedTemporaryFile() as existing_file:
            required_options = [
                existing_file.name,
                "--output-h5ad-file",
                existing_file.name,
                "--overwrite",
                "--model-url",
                model_file_path,
            ]
            result = CliRunner().invoke(
                annotate,
                required_options + [],
            )

            print(result.stdout)
            self.assertNotEqual(1, result.exit_code, "aborts with non-success code")
            self.assertIn(
                f"Wrote annotations to {existing_file.name}",
                result.stdout,
                "success message is correct on output file overwrite",
            )


# TODO:
# Test annotate cli args more comprehensively
# Test server.cli.annotate._validate_options
# Test model caching feature works
# Test model loading from s3 works (maybe w/just a real model)

if __name__ == "__main__":
    unittest.main()
