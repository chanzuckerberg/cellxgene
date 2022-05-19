import os.path
import unittest

import scanpy
from click.testing import CliRunner

from server.cli.annotate import annotate
from test.unit.annotate.fixtures.cell_type_annotate_model_fixture import build_fake_model, write_model, \
    write_query_dataset


class TestCliAnnotate(unittest.TestCase):

    def test__output_h5ad_file_option__writes_new_file(self):
        query_dataset_file_path = write_query_dataset(100, 10)
        model_file_path = write_model(build_fake_model({"x", "y", "z"}))

        result = CliRunner().invoke(annotate,
                                    ['--input-h5ad-file', query_dataset_file_path,
                                     '--model-url', model_file_path,
                                     '--output-h5ad-file', f"{query_dataset_file_path}.output"])

        self.assertEqual(0, result.exit_code, "runs successfully")
        self.assertTrue(os.path.isfile(f"{query_dataset_file_path}.output"), "output file written")

    def test__update_h5ad_file_option__updates_input_file(self):
        query_dataset_file_path = write_query_dataset(100, 10)
        model_file_path = write_model(build_fake_model({"x", "y", "z"}))

        result = CliRunner().invoke(annotate,
                                    ['--input-h5ad-file', query_dataset_file_path,
                                     '--model-url', model_file_path,
                                     '--update-h5ad-file'])

        self.assertEqual(0, result.exit_code, "runs successfully")
        self.assertIsNotNone(scanpy.read_h5ad(query_dataset_file_path).obs['cxg_predicted_cell_type'])


if __name__ == '__main__':
    unittest.main()
