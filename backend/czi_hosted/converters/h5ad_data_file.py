import json
import logging
from os import path

import anndata
import numpy as np
import tiledb

from backend.common_utils.colors import convert_anndata_category_colors_to_cxg_category_colors
from backend.czi_hosted.common.corpora import corpora_get_props_from_anndata
from backend.common_utils.errors import ColorFormatException
from backend.czi_hosted.common.utils.cxg_constants import CxgConstants
from backend.czi_hosted.common.utils.cxg_generation_utils import (
    convert_dictionary_to_cxg_group,
    convert_dataframe_to_cxg_array,
    convert_ndarray_to_cxg_dense_array,
    convert_matrix_to_cxg_array,
)
from backend.czi_hosted.common.utils.matrix_utils import is_matrix_sparse, get_column_shift_encode_for_matrix


class H5ADDataFile:
    """ Class encapsulating required information about an H5AD datafile that ultimately will be transformed into
    another format (currently just CXG is supported). """

    def __init__(
        self,
        input_filename,
        backed=False,
        dataset_title=None,
        dataset_about=None,
        obs_index_column_name=None,
        vars_index_column_name=None,
        use_corpora_schema=True,
    ):
        self.input_filename = input_filename
        self.backed = backed
        self.dataset_title = dataset_title
        self.dataset_about = dataset_about
        self.obs_index_column_name = obs_index_column_name
        self.vars_index_column_name = vars_index_column_name

        self.use_corpora_schema = use_corpora_schema

        self.validate_input_file_type()

        self.extract_anndata_elements_from_file()
        self.extract_metadata_about_dataset()

        self.validate_anndata()

    def to_cxg(self, output_cxg_directory, sparse_threshold, convert_anndata_colors_to_cxg_colors=True):
        """
        Writes the following attributes of the anndata to CXG: 1) the metadata as metadata attached to an empty
        DenseArray, 2) the obs DataFrame as a DenseArray, 3) the var DataFrame as a DenseArray, 4) all valid
        embeddings stored in obsm, each one as a DenseArray, 5) the main X matrix of the anndata as either a
        SparseArray or DenseArray based on the `sparse_threshold`, and optionally 6) the column shift of the main X
        matrix that might turn an otherwise Dense matrix into a Sparse matrix.
        """

        logging.info("Beginning writing to CXG.")
        ctx = tiledb.Ctx(
            {
                "sm.num_reader_threads": 32,
                "sm.num_writer_threads": 32,
                "sm.consolidation.buffer_size": 1 * 1024 * 1024 * 1024,
            }
        )

        tiledb.group_create(output_cxg_directory, ctx=ctx)
        logging.info(f"\t...group created, with name {output_cxg_directory}")

        convert_dictionary_to_cxg_group(
            output_cxg_directory, self.generate_cxg_metadata(convert_anndata_colors_to_cxg_colors)
        )
        logging.info("\t...dataset metadata saved")

        convert_dataframe_to_cxg_array(output_cxg_directory, "obs", self.obs, self.obs_index_column_name, ctx)
        logging.info("\t...dataset obs dataframe saved")

        convert_dataframe_to_cxg_array(output_cxg_directory, "var", self.var, self.var_index_column_name, ctx)
        logging.info("\t...dataset var dataframe saved")

        self.write_anndata_embeddings_to_cxg(output_cxg_directory, ctx)
        logging.info("\t...dataset embeddings saved")

        self.write_anndata_x_matrix_to_cxg(output_cxg_directory, ctx, sparse_threshold)
        logging.info("\t...dataset X matrix saved")

        logging.info("Completed writing to CXG.")

    def write_anndata_x_matrix_to_cxg(self, output_cxg_directory, ctx, sparse_threshold):
        matrix_container = f"{output_cxg_directory}/X"

        x_matrix_data = self.anndata.X
        is_sparse = is_matrix_sparse(x_matrix_data, sparse_threshold)
        if not is_sparse:
            col_shift = get_column_shift_encode_for_matrix(x_matrix_data, sparse_threshold)
            is_sparse = col_shift is not None
        else:
            col_shift = None

        if col_shift is not None:
            logging.info("Converting matrix X as sparse matrix with column shift encoding")
            x_col_shift_name = f"{output_cxg_directory}/X_col_shift"
            convert_ndarray_to_cxg_dense_array(x_col_shift_name, col_shift, ctx)

        convert_matrix_to_cxg_array(matrix_container, x_matrix_data, is_sparse, ctx, col_shift)

        tiledb.consolidate(matrix_container, ctx=ctx)
        if hasattr(tiledb, "vacuum"):
            tiledb.vacuum(matrix_container)

    def write_anndata_embeddings_to_cxg(self, output_cxg_directory, ctx):
        def is_valid_embedding(adata, embedding_name, embedding_array):
            """
            Returns true if this layout data is a valid array for front-end presentation with the following criteria:
                * ndarray, with shape (n_obs, >= 2), dtype float/int/uint
                * follows ScanPy embedding naming conventions
                * with all values finite or NaN (no +Inf or -Inf)
            """

            is_valid = isinstance(embedding_name, str) and embedding_name.startswith("X_") and len(embedding_name) > 2
            is_valid = is_valid and isinstance(embedding_array, np.ndarray) and embedding_array.dtype.kind in "fiu"
            is_valid = is_valid and embedding_array.shape[0] == adata.n_obs and embedding_array.shape[1] >= 2
            is_valid = is_valid and not np.any(np.isinf(embedding_array)) and not np.all(np.isnan(embedding_array))
            return is_valid

        embedding_container = f"{output_cxg_directory}/emb"
        tiledb.group_create(embedding_container, ctx=ctx)

        for embedding_name, embedding_values in self.anndata.obsm.items():
            if is_valid_embedding(self.anndata, embedding_name, embedding_values):
                embedding_name = f"{embedding_container}/{embedding_name[2:]}"
                convert_ndarray_to_cxg_dense_array(embedding_name, embedding_values, ctx)
                logging.info(f"\t\t...{embedding_name} embedding created")

    def generate_cxg_metadata(self, convert_anndata_colors_to_cxg_colors):
        """
        Return a dictionary containing metadata about CXG dataset. This include data about the version as well as
        Corpora schema properties if they exist, among other pieces of metadata.
        """

        cxg_group_metadata = {
            "cxg_version": CxgConstants.CXG_VERSION,
            "cxg_properties": json.dumps({"title": self.dataset_title, "about": self.dataset_about}),
        }
        if self.corpora_properties is not None:
            cxg_group_metadata["corpora"] = json.dumps(self.corpora_properties)

        if convert_anndata_colors_to_cxg_colors:
            try:
                cxg_group_metadata["cxg_category_colors"] = json.dumps(
                    convert_anndata_category_colors_to_cxg_category_colors(self.anndata)
                )
            except ColorFormatException:
                logging.warning(
                    "Failed to extract colors from H5AD file! Fix the H5AD file or rerun with "
                    "--disable-custom-colors. See help for more details."
                )

        return cxg_group_metadata

    def validate_input_file_type(self):
        """
        Validate that the input file is of a type that we can handle. Currently the only valid file type is `.h5ad`.
        """

        if not self.input_filename.endswith(".h5ad"):
            raise Exception(f"Cannot process input file {self.input_filename}. File must be an H5AD.")

        if self.dataset_title or self.dataset_about:
            logging.warning(
                "If you convert this dataset into CXG and you explicit specify values for the dataset title metadata "
                "or the dataset about metadata, it will override any metadata that is extracted as part of the "
                "Corpora schema fields."
            )

    def validate_anndata(self):
        if not self.var.index.is_unique:
            raise ValueError("Variable index in AnnData object is not unique.")
        if not self.obs.index.is_unique:
            raise ValueError("Observation index in AnnData object is not unique.")

    def extract_anndata_elements_from_file(self):
        logging.info(f"Reading in AnnData dataset: {path.basename(self.input_filename)}")
        self.anndata = anndata.read_h5ad(self.input_filename, backed="r" if self.backed else None)
        logging.info("Completed reading in AnnData dataset!")

        self.obs = self.transform_dataframe_index_into_column(self.anndata.obs, "obs", self.obs_index_column_name)
        self.var = self.transform_dataframe_index_into_column(self.anndata.var, "var", self.vars_index_column_name)

    def extract_metadata_about_dataset(self):
        """
        Extract metadata information about the dataset that upon conversion will be saved as group metadata with the
        CXG that is generated. This metadata information includes Corpora schema properties, the dataset title and
        a link that details more information about the dataset.
        """

        self.corpora_properties = corpora_get_props_from_anndata(self.anndata) if self.use_corpora_schema else None
        if self.corpora_properties is None and self.use_corpora_schema:
            # If the return value is None, this means that we were not able to figure out what version of the Corpora
            # schema the object is using and therefore cannot extract any properties.
            raise ValueError("Unknown source file schema version is unsupported.")

        # The title and about properties of the dataset are set by the following order: if they are explicitly defined
        # then use the explicit value. If the dataset is a Corpora-schema based schema, then extract the title and about
        # from the corpora_properties. Otherwise, use the input filename (only for title, about will be blank).
        if self.corpora_properties:
            corpora_project_links = self.corpora_properties.get("project_links", [])
            corpora_about_link = next(
                (link for link in corpora_project_links if (link.get("link_type", None) == "SUMMARY")), {}
            )
        else:
            corpora_about_link = {}

        filename = path.splitext(path.basename(self.input_filename))[0]

        self.dataset_title = self.dataset_title if self.dataset_title else corpora_about_link.get("link_name", filename)
        self.dataset_about = self.dataset_about if self.dataset_about else corpora_about_link.get("link_url")

    def transform_dataframe_index_into_column(self, dataframe, dataframe_name, index_column_name):
        """
        Convert the dataframe's index into another column in the dataframe. If an index_column_name is specified,
        use that column as the index instead.
        """

        if index_column_name is None:
            # Create a unique column name for the index.
            suffix = 0
            while f"name_{suffix}" in dataframe.columns:
                suffix += 1
            index_column_name = f"name_{suffix}"

            # Turn the index into a normal column
            dataframe.rename_axis(index_column_name, inplace=True)
            dataframe.reset_index(inplace=True)

        elif index_column_name in dataframe.columns:
            # User has specified alternative column for unique names, and it exists
            if not dataframe[index_column_name].is_unique:
                raise KeyError(
                    f"Values in {dataframe_name}.{index_column_name} must be unique. Please prepare data to contain "
                    f"unique values."
                )
        else:
            raise KeyError(f"Column {index_column_name} does not exist.")

        setattr(self, f"{dataframe_name}_index_column_name", index_column_name)
        return dataframe
