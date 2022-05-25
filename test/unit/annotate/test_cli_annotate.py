import os.path
import unittest

import scanpy
from click.testing import CliRunner

from server.cli.annotate import annotate
from test.unit.annotate.fixtures.cell_type_annotate_model_fixture import write_model, \
    write_query_dataset, FakeModel, build_dataset


class TestCliAnnotate(unittest.TestCase):

    def test__local_model_file__loads(self):
        query_dataset_file_path = write_query_dataset(100, gene_identifiers=['a', 'b', 'c'])
        labels = {"x", "y", "z"}
        ref_dataset = build_dataset(1000, gene_identifiers=['a', 'b', 'c'])
        model_file_path = write_model(FakeModel(labels, ref_dataset))

        result = CliRunner().invoke(annotate,
                                    ['--input-h5ad-file', query_dataset_file_path,
                                     '--model-url', model_file_path,
                                     '--output-h5ad-file', f"{query_dataset_file_path}.output"])

        self.assertEqual(0, result.exit_code, "runs successfully")

    @unittest.skipUnless(os.getenv('TEST_REMOTE_MODEL_URL'), 'remote model url env var missing')
    # For testing with a real model stores in a remote location (e.g. S3), if available. We don't want to setup a fake
    # S3 service, so this test is intended to be run in an ad hoc fashion if related code changes have been made.
    def test__s3_model_file__loads(self):
        query_dataset_file_path = write_query_dataset(100, 10)

        result = CliRunner().invoke(annotate,
                                    ['--input-h5ad-file', query_dataset_file_path,
                                     '--model-url', os.getenv('TEST_REMOTE_MODEL_URL'),
                                     '--output-h5ad-file', f"{query_dataset_file_path}.output"])

        self.assertEqual(0, result.exit_code, "runs successfully")

    def test__update_h5ad_file_option__updates_input_file(self):
        query_dataset_file_path = write_query_dataset(100, 10)
        labels = {"x", "y", "z"}
        model_file_path = write_model(FakeModel(labels))

        result = CliRunner().invoke(annotate,
                                    ['--input-h5ad-file', query_dataset_file_path,
                                     '--model-url', model_file_path,
                                     '--update-h5ad-file'])

        self.assertEqual(0, result.exit_code, "runs successfully")
        self.assertIsNotNone(scanpy.read_h5ad(query_dataset_file_path).obs['cxg_predicted_cell_type'])

    def test__annotation_column_options__writes_correct_output_column_name(self):
        query_dataset_file_path = write_query_dataset(100, 10)
        labels = {"x", "y", "z"}
        model_file_path = write_model(FakeModel(labels))

        result = CliRunner().invoke(annotate,
                                    ['--input-h5ad-file', query_dataset_file_path,
                                     '--model-url', model_file_path,
                                     '--update-h5ad-file',
                                     '--annotation-column-prefix', 'test_prefix',
                                     '--run-name', 'test_run'])

        self.assertEqual(0, result.exit_code, "runs successfully")
        self.assertIn('test_prefix_cell_type_test_run', scanpy.read_h5ad(query_dataset_file_path).obs.keys())

    # TODO: Test S3 model locations (probably not possible)


if __name__ == '__main__':
    unittest.main()
