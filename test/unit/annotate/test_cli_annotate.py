import os.path
import unittest
from unittest import mock
from unittest.mock import patch

import scanpy
from click.testing import CliRunner

from server.cli.annotate import annotate
from test.unit.annotate.fixtures.cell_type_annotate_model_fixture import write_model, \
    write_query_dataset, FakeModel, build_dataset


class TestCliAnnotate(unittest.TestCase):

    def test__local_model_file__loads(self):
        query_dataset_file_path = write_query_dataset(100, gene_identifiers=['a', 'b', 'c'])
        ref_dataset = build_dataset(1000, gene_identifiers=['a', 'b', 'c'])
        model_file_path = write_model(FakeModel(ref_dataset=ref_dataset))

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
        query_dataset_file_path = write_query_dataset(100)
        model_file_path = write_model(FakeModel())

        result = CliRunner().invoke(annotate,
                                    ['--input-h5ad-file', query_dataset_file_path,
                                     '--model-url', model_file_path,
                                     '--update-h5ad-file'])

        print(result.output)
        self.assertEqual(0, result.exit_code, "runs successfully")
        self.assertIsNotNone(scanpy.read_h5ad(query_dataset_file_path).obs['cxg_predicted_cell_type'])

    def test__annotation_column_options__writes_correct_output_column_name(self):
        query_dataset_file_path = write_query_dataset(100)
        model_file_path = write_model(FakeModel())

        result = CliRunner().invoke(annotate,
                                    ['--input-h5ad-file', query_dataset_file_path,
                                     '--model-url', model_file_path,
                                     '--update-h5ad-file',
                                     '--annotation-column-prefix', 'test_prefix',
                                     '--run-name', 'test_run'])

        self.assertEqual(0, result.exit_code, "runs successfully")
        self.assertIn('test_prefix_cell_type_test_run', scanpy.read_h5ad(query_dataset_file_path).obs.keys())

    def test__min_common_gene_pct_option__aborts_if_violated(self):
        query_dataset_file_path = write_query_dataset(10, gene_identifiers=[f'gene_{i}' for i in range(49)])
        ref_dataset = build_dataset(10, gene_identifiers=[f'gene_{i}' for i in range(100)])
        model_file_path = write_model(FakeModel(ref_dataset=ref_dataset))

        result = CliRunner().invoke(annotate,
                                    ['--input-h5ad-file', query_dataset_file_path,
                                     '--model-url', model_file_path,
                                     '--update-h5ad-file',
                                     '--min-common-gene-pct', '50',
                                     '--run-name', 'test_run'])

        self.assertNotEqual(0, result.exit_code, "aborts if --min-common-gene-pct violated")
        self.assertEqual(
            'There are only 49 genes in common between the query and reference datasets, but at least 50 are required',
            str(result.exception))


    def test__min_common_gene_pct_option__succeeds_if_boundary_value_satisfied(self):
        query_dataset_file_path = write_query_dataset(10, gene_identifiers=[f'gene_{i}' for i in range(51)])
        ref_dataset = build_dataset(10, gene_identifiers=[f'gene_{i}' for i in range(100)])
        model_file_path = write_model(FakeModel(ref_dataset=ref_dataset))

        result = CliRunner().invoke(annotate,
                                    ['--input-h5ad-file', query_dataset_file_path,
                                     '--model-url', model_file_path,
                                     '--update-h5ad-file',
                                     '--min-common-gene-pct', '50',
                                     '--run-name', 'test_run'])

        self.assertEqual(0, result.exit_code, "succeeds if --min-common-gene-pct satisfied")

    def test__min_common_gene_pct_option__validates_value(self):
        required_options = ['--input-h5ad-file', 'dummy',
                            '--model-url', 'dummy',
                            '--update-h5ad-file',
                            '--run-name', 'test_run']

        result = CliRunner().invoke(annotate, required_options + ['--min-common-gene-pct', '0'])

        self.assertNotEqual(0, result.exit_code, "validates --min-common-gene-pct > 0")
        self.assertIn("Error: Invalid value for '--min-common-gene-pct': 0 is not in the range 0<x<=100.",
                      result.output)

        result = CliRunner().invoke(annotate, required_options + ['--min-common-gene-pct', '101'])

        self.assertNotEqual(0, result.exit_code, "validates --min-common-gene-pct <= 100")
        self.assertIn("Error: Invalid value for '--min-common-gene-pct': 101 is not in the range 0<x<=100.",
                      result.output)

    def test__adds_metadata(self):
        query_dataset_file_path = write_query_dataset(1)

        mock_fetch = mock.Mock()
        mock_fetch.return_value = FakeModel()
        with patch('server.cli.annotate._fetch_model', mock_fetch):
            result = CliRunner().invoke(annotate,
                                        ['--input-h5ad-file', str(query_dataset_file_path),
                                         '--model-url', 'http://model_respository.org/annotation/cell_type/model1.pkl',
                                         '--update-h5ad-file',
                                         '--run-name', 'test_run'])

        self.assertEqual(0, result.exit_code)  # sanity check

        ad = scanpy.read_h5ad(query_dataset_file_path)
        self.assertEqual('http://model_respository.org/annotation/cell_type/model1.pkl',
                         ad.uns['cxg_predictions']['cxg_predicted_cell_type_test_run']['model_url'])

    def test__metadata__appends_existing(self):
        query_dataset_file_path = write_query_dataset(1)
        ad = scanpy.read_h5ad(str(query_dataset_file_path))
        ad.uns['existing_metadata'] = {"keep_me": "please"}
        ad.write_h5ad(query_dataset_file_path)

        result = CliRunner().invoke(annotate,
                                    ['--input-h5ad-file', str(query_dataset_file_path),
                                     '--model-url', str(write_model(FakeModel())),
                                     '--update-h5ad-file',
                                     '--run-name', 'test_run'])

        self.assertEqual(0, result.exit_code)  # sanity check

        ad = scanpy.read_h5ad(query_dataset_file_path)
        self.assertTrue('existing_metadata' in ad.uns)


if __name__ == '__main__':
    unittest.main()
