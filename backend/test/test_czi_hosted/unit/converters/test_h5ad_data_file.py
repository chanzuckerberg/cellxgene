import json
import unittest
from glob import glob
from os import remove, path
from shutil import rmtree
from uuid import uuid4

import anndata
import numpy as np
from pandas import Series, DataFrame
import tiledb
from urllib.parse import urljoin

from backend.czi_hosted.common.corpora import CorporaConstants
from backend.czi_hosted.converters.h5ad_data_file import H5ADDataFile

from backend.test import PROJECT_ROOT


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

        self._validate_cxg_and_h5ad_content_match(self.sample_h5ad_filename, self.sample_output_directory, True)

    def test__to_cxg__simple_anndata_with_corpora_and_sparse(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename)
        h5ad_file.to_cxg(self.sample_output_directory, 100)

        self._validate_cxg_and_h5ad_content_match(self.sample_h5ad_filename, self.sample_output_directory, True)

    def test__to_cxg__simple_anndata_no_corpora_and_dense(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename, use_corpora_schema=False)
        h5ad_file.to_cxg(self.sample_output_directory, 0)

        self._validate_cxg_and_h5ad_content_match(self.sample_h5ad_filename, self.sample_output_directory, False)

    def test__to_cxg__simple_anndata_with_corpora_and_dense(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename)
        h5ad_file.to_cxg(self.sample_output_directory, 0)

        self._validate_cxg_and_h5ad_content_match(self.sample_h5ad_filename, self.sample_output_directory, False)

    def test__to_cxg__with_sparse_column_encoding(self):
        anndata = self._create_sample_anndata_dataset()
        anndata.X = np.ones((3, 4))
        sparse_with_column_shift_filename = self._write_anndata_to_file(anndata)

        h5ad_file = H5ADDataFile(sparse_with_column_shift_filename)
        h5ad_file.to_cxg(self.sample_output_directory, 50)

        self._validate_cxg_and_h5ad_content_match(
            sparse_with_column_shift_filename, self.sample_output_directory, False, has_column_encoding=True
        )

        # Clean up
        remove(sparse_with_column_shift_filename)

    def _validate_cxg_and_h5ad_content_match(self, h5ad_filename, cxg_directory, is_sparse, has_column_encoding=False):
        anndata_object = anndata.read_h5ad(h5ad_filename)
        tiledb_ctx = tiledb.Ctx(
            {"sm.tile_cache_size": 8 * 1024 * 1024 * 1024, "sm.num_reader_threads": 32, "vfs.s3.region": "us-east-1"}
        )

        # Array locations
        metadata_array_location = f"{cxg_directory}/cxg_group_metadata"
        main_x_array_location = f"{cxg_directory}/X"
        embedding_array_location = f"{cxg_directory}/emb"
        specific_embedding_array_location = f"{self.sample_output_directory}/emb/awesome_embedding"
        obs_array_location = f"{cxg_directory}/obs"
        var_array_location = f"{cxg_directory}/var"
        x_col_shift_array_location = f"{cxg_directory}/X_col_shift"

        # Assert CXG structure
        self.assertEqual(tiledb.object_type(cxg_directory, ctx=tiledb_ctx), "group")
        self.assertEqual(tiledb.object_type(obs_array_location, ctx=tiledb_ctx), "array")
        self.assertEqual(tiledb.object_type(var_array_location, ctx=tiledb_ctx), "array")
        self.assertEqual(tiledb.object_type(main_x_array_location, ctx=tiledb_ctx), "array")
        self.assertEqual(tiledb.object_type(embedding_array_location, ctx=tiledb_ctx), "group")
        self.assertEqual(tiledb.object_type(specific_embedding_array_location, ctx=tiledb_ctx), "array")

        if has_column_encoding:
            self.assertEqual(tiledb.object_type(x_col_shift_array_location, ctx=tiledb_ctx), "array")

        # Validate metadata
        metadata_array = tiledb.DenseArray(metadata_array_location, mode="r", ctx=tiledb_ctx)
        self.assertIn("cxg_version", metadata_array.meta)

        # Validate obs index
        obs_array = tiledb.DenseArray(obs_array_location, mode="r", ctx=tiledb_ctx)
        expected_index_data = anndata_object.obs.index.to_numpy()
        index_name = json.loads(obs_array.meta["cxg_schema"])["index"]
        actual_index_data = obs_array.query(attrs=[index_name])[:][index_name]
        self.assertTrue(np.array_equal(expected_index_data, actual_index_data))

        # Validate obs columns
        expected_columns = list(anndata_object.obs.columns.values)
        for column_name in expected_columns:
            expected_data = anndata_object.obs[column_name].to_numpy()
            actual_data = obs_array.query(attrs=[column_name])[:][column_name]
            self.assertTrue(np.array_equal(expected_data, actual_data))

        # Validate var index
        var_array = tiledb.DenseArray(var_array_location, mode="r", ctx=tiledb_ctx)
        expected_index_data = anndata_object.var.index.to_numpy()
        index_name = json.loads(var_array.meta["cxg_schema"])["index"]
        actual_index_data = var_array.query(attrs=[index_name])[:][index_name]
        self.assertTrue(np.array_equal(expected_index_data, actual_index_data))

        # Validate var columns
        expected_columns = list(anndata_object.var.columns.values)
        for column_name in expected_columns:
            expected_data = anndata_object.var[column_name].to_numpy()
            actual_data = var_array.query(attrs=[column_name])[:][column_name]
            self.assertTrue(np.array_equal(expected_data, actual_data))

        # Validate embedding
        expected_embedding_data = anndata_object.obsm.get("X_awesome_embedding")
        embedding_array = tiledb.DenseArray(specific_embedding_array_location, mode="r", ctx=tiledb_ctx)
        actual_embedding_data = embedding_array[:, 0:2]
        self.assertTrue(np.array_equal(expected_embedding_data, actual_embedding_data))

        # Validate X matrix if not column shifted
        if not has_column_encoding:
            expected_x_data = anndata_object.X
            if is_sparse:
                x_array = tiledb.SparseArray(main_x_array_location, mode="r", ctx=tiledb_ctx)
                actual_x_data = np.reshape(x_array[:, :][""], expected_x_data.shape)
            else:
                x_array = tiledb.DenseArray(main_x_array_location, mode="r", ctx=tiledb_ctx)
                actual_x_data = x_array[:, :]
            self.assertTrue(np.array_equal(expected_x_data, actual_x_data))

    def _write_anndata_to_file(self, anndata):
        temporary_filename = f"{PROJECT_ROOT}/backend/test/fixtures/{uuid4()}.h5ad"
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
