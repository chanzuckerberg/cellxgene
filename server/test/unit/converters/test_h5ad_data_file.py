import json
import unittest
from glob import glob
from os import popen, remove, path
from shutil import rmtree
from uuid import uuid4

import anndata
import numpy as np
from pandas import Series, DataFrame

from server.common.utils.corpora_constants import CorporaConstants
from server.converters.h5ad_data_file import H5ADDataFile

PROJECT_ROOT = popen("git rev-parse --show-toplevel").read().strip()


class TestH5ADDataFile(unittest.TestCase):
    def setUp(self):
        self.sample_anndata = self._create_sample_anndata_dataset()
        self.sample_h5ad_filename = self._write_anndata_to_file(self.sample_anndata)

        self.sample_output_directory = path.splitext(self.sample_h5ad_filename)[0] + ".cxg"

    def tearDown(self):
        if self.sample_h5ad_filename:
            remove(self.sample_h5ad_filename)

        if path.isdir(self.sample_output_directory):
            rmtree(self.sample_output_directory)

    def test__create_h5ad_data_file__non_h5ad_raises_exception(self):
        non_h5ad_filename = "my_fancy_dataset.csv"

        with self.assertRaises(Exception) as exception_context:
            H5ADDataFile(non_h5ad_filename)

        self.assertIn("File must be an H5AD", str(exception_context.exception))

    def test__create_h5ad_data_file__assert_warning_outputted_if_dataset_title_or_about_given(self):
        with self.assertLogs(level="WARN") as logger:
            H5ADDataFile(
                self.sample_h5ad_filename,
                dataset_title="My Awesome Dataset",
                dataset_about="http://www.awesomedataset.com",
                use_corpora_schema=False,
            )

        self.assertIn("will override any metadata that is extracted", logger.output[0])

    def test__create_h5ad_data_file__reads_anndata_successfully(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename, use_corpora_schema=False)

        self.assertTrue((h5ad_file.anndata.X == self.sample_anndata.X).all())
        self.assertEqual(
            h5ad_file.anndata.obs.sort_index(inplace=True), self.sample_anndata.obs.sort_index(inplace=True)
        )
        self.assertEqual(
            h5ad_file.anndata.var.sort_index(inplace=True), self.sample_anndata.var.sort_index(inplace=True)
        )

        for key in h5ad_file.anndata.obsm.keys():
            self.assertIn(key, self.sample_anndata.obsm.keys())
            self.assertTrue((h5ad_file.anndata.obsm[key] == self.sample_anndata.obsm[key]).all())

        for key in self.sample_anndata.obsm.keys():
            self.assertIn(key, h5ad_file.anndata.obsm.keys())
            self.assertTrue((h5ad_file.anndata.obsm[key] == self.sample_anndata.obsm[key]).all())

    def test__create_h5ad_data_file__copies_index_of_obs_and_var_to_column(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename, use_corpora_schema=False)

        # The automatic name chosen for the index should be "name_0"
        self.assertNotIn("name_0", self.sample_anndata.obs.columns)
        self.assertIn("name_0", h5ad_file.obs.columns)

        self.assertNotIn("name_0", self.sample_anndata.var.columns)
        self.assertIn("name_0", h5ad_file.var.columns)

    def test__create_h5ad_data_file__no_copy_if_obs_and_var_index_names_specified(self):
        h5ad_file = H5ADDataFile(
            self.sample_h5ad_filename,
            use_corpora_schema=False,
            obs_index_column_name="float_category",
            vars_index_column_name="int_category",
        )

        self.assertNotIn("name_0", h5ad_file.obs.columns)
        self.assertNotIn("name_0", h5ad_file.var.columns)

    def test__create_h5ad_data_file__obs_and_var_index_names_specified_not_unique_raises_exception(self):

        with self.assertRaises(Exception) as exception_context:
            H5ADDataFile(
                self.sample_h5ad_filename,
                use_corpora_schema=False,
                obs_index_column_name="float_category",
                vars_index_column_name="bool_category",
            )

        self.assertIn("Please prepare data to contain unique values", str(exception_context.exception))

    def test__create_h5ad_data_file__obs_and_var_index_names_specified_doesnt_exist_raises_exception(self):
        with self.assertRaises(Exception) as exception_context:
            H5ADDataFile(
                self.sample_h5ad_filename,
                use_corpora_schema=False,
                obs_index_column_name="unknown_category",
                vars_index_column_name="i_dont_exist",
            )

        self.assertIn("does not exist", str(exception_context.exception))

    def test__create_h5ad_data_file__extract_about_and_title_from_dataset(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename)

        self.assertEqual(h5ad_file.dataset_title, "random_link_name")
        self.assertEqual(h5ad_file.dataset_about, "www.link.com")

    def test__create_h5ad_data_file__inputted_dataset_title_and_about_overrides_extracted(self):
        h5ad_file = H5ADDataFile(
            self.sample_h5ad_filename, dataset_about="override_about", dataset_title="override_title"
        )

        self.assertEqual(h5ad_file.dataset_title, "override_title")
        self.assertEqual(h5ad_file.dataset_about, "override_about")

    def test__to_cxg__simple_anndata_no_corpora_and_sparse(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename, use_corpora_schema=False)
        h5ad_file.to_cxg(self.sample_output_directory, 100)

        self._validate_expected_generated_list_of_tiledb_files()

    def test__to_cxg__simple_anndata_with_corpora_and_sparse(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename)
        h5ad_file.to_cxg(self.sample_output_directory, 100)

        self._validate_expected_generated_list_of_tiledb_files()

    def test__to_cxg__simple_anndata_no_corpora_and_dense(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename, use_corpora_schema=False)
        h5ad_file.to_cxg(self.sample_output_directory, 0)

        self._validate_expected_generated_list_of_tiledb_files()

    def test__to_cxg__simple_anndata_with_corpora_and_dense(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename)
        h5ad_file.to_cxg(self.sample_output_directory, 0)

        self._validate_expected_generated_list_of_tiledb_files()

    def test__to_cxg__with_sparse_column_encoding(self):
        anndata = self._create_sample_anndata_dataset()
        anndata.X = np.ones((3, 4))
        sparse_with_column_shift_filename = self._write_anndata_to_file(anndata)

        h5ad_file = H5ADDataFile(sparse_with_column_shift_filename)
        h5ad_file.to_cxg(self.sample_output_directory, 50)

        self._validate_expected_generated_list_of_tiledb_files(has_column_encoding=True)

        # Clean up
        remove(sparse_with_column_shift_filename)

    def _validate_expected_generated_list_of_tiledb_files(self, has_column_encoding=False):
        (
            expected_directories,
            expected_obs_files,
            expected_var_files,
        ) = self._get_expected_generated_list_of_tiledb_files()

        for directory in expected_directories:
            self.assertTrue(path.isdir(directory))

        for obs_file in expected_obs_files:
            expected_location_of_obs_file = f"{self.sample_output_directory}/obs/*/{obs_file}"
            self.assertTrue(path.isfile(glob(expected_location_of_obs_file)[0]))

        for var_file in expected_var_files:
            expected_location_of_var_file = f"{self.sample_output_directory}/var/*/{var_file}"
            self.assertTrue(path.isfile(glob(expected_location_of_var_file)[0]))

        if has_column_encoding:
            self.assertTrue(path.isdir(f"{self.sample_output_directory}/X_col_shift"))

    def _get_expected_generated_list_of_tiledb_files(self):

        # Expected directories
        metadata_directory = f"{self.sample_output_directory}/cxg_group_metadata"
        main_x_directory = f"{self.sample_output_directory}/X"
        overall_embedding_directory = f"{self.sample_output_directory}/emb"
        specific_embedding_directory = f"{self.sample_output_directory}/emb/awesome_embedding"
        obs_directory = f"{self.sample_output_directory}/obs"
        var_directory = f"{self.sample_output_directory}/var"

        # Obs files
        obs_files = []
        obs_files.append("name_0.tdb")
        obs_files.append("name_0_var.tdb")
        obs_files.append("string_category.tdb")
        obs_files.append("string_category_var.tdb")
        obs_files.append("float_category.tdb")

        # Var files
        var_files = []
        var_files.append("name_0.tdb")
        var_files.append("name_0_var.tdb")
        var_files.append("bool_category.tdb")
        var_files.append("int_category.tdb")

        return (
            [
                metadata_directory,
                main_x_directory,
                overall_embedding_directory,
                specific_embedding_directory,
                obs_directory,
                var_directory,
            ],
            obs_files,
            var_files,
        )

    def _write_anndata_to_file(self, anndata):
        temporary_filename = f"{PROJECT_ROOT}/server/test/fixtures/{uuid4()}.h5ad"
        anndata.write(temporary_filename)

        return temporary_filename

    def _create_sample_anndata_dataset(self):
        # Create X
        X = np.random.rand(3, 4)

        # Create obs
        random_string_category = Series(data=["a", "b", "b"], dtype="category")
        random_float_category = Series(data=[3.2, 1.1, 2.2], dtype=np.float32)
        obs_dataframe = DataFrame(
            data={"string_category": random_string_category, "float_category": random_float_category}
        )
        obs = obs_dataframe

        # Create vars
        random_int_category = Series(data=[3, 1, 2, 4], dtype=np.int32)
        random_bool_category = Series(data=[True, True, False, True], dtype=np.bool_)
        var_dataframe = DataFrame(data={"int_category": random_int_category, "bool_category": random_bool_category})
        var = var_dataframe

        # Create embeddings
        random_embedding = np.random.rand(3, 2)
        obsm = {"X_awesome_embedding": random_embedding}

        # Create uns corpora metadata
        uns = {}
        for metadata_field in CorporaConstants.REQUIRED_SIMPLE_METADATA_FIELDS:
            uns[metadata_field] = "random"

        for metadata_field in CorporaConstants.OPTIONAL_JSON_ENCODED_METADATA_FIELD:
            uns[metadata_field] = json.dumps({"random_key": "random_value"})

        # Need to carefully set the corpora schema versions in order for tests to pass.
        uns["version"] = {"corpora_schema_version": "1.0.0", "corpora_encoding_version": "0.1.0"}

        # Set project links to be a dictionary
        uns["project_links"] = json.dumps(
            [{"link_name": "random_link_name", "link_url": "www.link.com", "link_type": "SUMMARY"}]
        )

        return anndata.AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=uns)
