import unittest
from tempfile import mkstemp

from click.testing import CliRunner

from server.cli.annotate import annotate
from test.unit.cli.mlflow_model_fixture import FakeModel, write_model


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

        query_dataset_file_path = mkstemp()
        model = FakeModel()
        model_file_path = write_model(model)

        result = CliRunner().invoke(
            annotate,
            [
                "--input-h5ad-file",
                query_dataset_file_path,
                "--model-url",
                model_file_path,
                "--output-h5ad-file",
                f"{query_dataset_file_path}.output",
            ],
        )

        # to help debugging, show the output from the CliRunner and MLflow stdout
        if result.exit_code:
            print(result.stdout)

        self.assertEqual(0, result.exit_code, "runs successfully")

        # The FakeModel will print it inputs to stdout, as "MODEL_INPUT={...}", allowing us to assert that it received valid inputs.
        self.assertIn(
            "MODEL_INPUT={"
            f'"query_dataset_h5ad_path": "{query_dataset_file_path}", '
            f'"output_h5ad_path": "{query_dataset_file_path}.output", '
            '"annotation_prefix": "cxg_cell_type", "classifier": "fine", '
            '"organism": "Homo sapiens", "use_gpu": true}',
            result.stdout,
            "inputs passed correctly",
        )


# TODO:
# Test annotate cli args more comprehensively
# Test server.cli.annotate._validate_options
# Test model caching feature works
# Test model loading from s3 works (maybe w/just a real model)

if __name__ == "__main__":
    unittest.main()
