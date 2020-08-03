import json
import logging
import re
from os import path

import anndata
import click
import numpy as np
import tiledb
from scipy.stats import mode

from server.common.colors import convert_anndata_category_colors_to_cxg_category_colors
from server.common.corpora import corpora_get_props_from_anndata
from server.common.errors import ColorFormatException

# The CXG container version number.  Must be a semver string (major.minor.patch)
# DO NOT UPDATE THIS WITHOUT ALSO UPDATING THE CXG SPECIFICATION.
CXG_VERSION = "0.2.0"


@click.command(
    name="convert",
    short_help="Converts an H5AD dataset to the CXG format.",
    help="Converts an H5AD dataset to the CXG format. The CXG format is a cellxgene-private data format "
    "that has performance and access characteristics amenable to a multi-dataset, multi-user serving "
    "environment. You will be able to launch the cellxgene using the `cellxgene launch` command as "
    "usually with the generated CXG file.",
)
@click.argument(
    "input-file",
    nargs=1,
    help="Path to the H5AD input file to be converted.",
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "-o",
    "--output-dir",
    help="Name of the output CXG directory. If not provided, will default to be the input filename with a "
    "CXG extension.",
)
@click.option(
    "-b",
    "--backed",
    help="When true, loads the H5AD in file backed mode. This will cause the conversion to be slower, "
    "but will use less memory.",
    default=False,
    show_default=True,
    is_flag=True,
)
@click.option(
    "-t",
    "--title",
    help="Human readable dataset title that will be included as metadata about the CXG file. If omitted, "
    "the dataset title will be the filename.",
)
@click.option(
    "-a",
    "--about",
    help="A fully qualified URL that provides more information about the dataset and will be included as "
    "metadata about the CXG file.",
)
@click.option(
    "-s",
    "--sparse-threshold",
    help="If the dataset's percent of non-zero values falls belows the specified threshold, then the X "
    "array of the dataset will be sparse. Since the default value is 0.0, the default will be to "
    "convert to dense array.",
    default=0.0,
    show_default=True,
)
@click.option("--obs-names", help="The name of the observation annotations.")  # TODO expand this
@click.option("--var-names", help="The name of the variables annotations.")  # TODO expand this
@click.option(
    "--disable-custom-colors",
    help="When set, conversion process will not extract scanpy-compatible category colors from the H5AD " "file.",
    default=False,
    show_default=True,
    is_flag=True,
)
@click.option(
    "--disable-corpora-schema",
    "When set, conversion process will neither extract nor store Corpora schema information. See "
    "https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema.md for "
    "more information.",
    default=False,
    show_default=True,
    is_flag=True,
)
@click.option(
    "--overwrite",
    help="When set to true, will overwrite the output file if the output file already exists.",
    default=False,
    show_default=True,
    is_flag=True,
)
@click.option("-v", "--verbose", count=True)
@click.help_option("--help", "-h", help="Show this message and exit.")
def convert_to_cxg(
    input_file,
    output_directory,
    backed,
    title,
    about,
    sparse_threshold,
    obs_names,
    var_names,
    disable_custom_colors,
    disable_corpora_schema,
    should_overwrite,
):
    """
    Convert a dataset file into CXG.
    """

    # Validate the command parameters
    validate_input_file_type(input_file)

    # Read in dataset
    logging.info(f"Reading in AnnData dataset: {path.basename(input_file)}")
    dataset = anndata.read_h5ad(input_file, backed="r" if backed else None)
    logging.info("Completed reading in AnnData dataset!")

    # Get the directory that will hold all the CXG files
    cxg_output_container = get_output_directory(input_file, output_directory, should_overwrite)

    # Get Corpora schema properties from the anndata object
    logging.info("Extracting metadata about the dataset")
    corpora_properties = extract_corpora_properties(dataset, title, about) if not disable_corpora_schema else None

    # Get the title and about properties of the dataset
    title, about = extract_title_and_about_dataset(title, about, corpora_properties, input_file)

    logging.info("Creating CXG metadata information")
    cxg_group_metadata = create_cxg_group_metadata(dataset, title, about, corpora_properties, not disable_custom_colors)

    logging.info("Beginning CXG dataset generation")
    generate_cxg_dataset(dataset, cxg_output_container, cxg_group_metadata, obs_names, var_names, sparse_threshold)
    logging.info("Completed CXG dataset generation!")


def generate_cxg_dataset(anndata_dataset, cxg_container, cxg_group_metadata, obs_names, var_names, sparse_threshold):
    # Clean and validate the AnnData dataset
    logging.info("Cleaning and validating AnnData object.")
    validate_anndata_object(anndata_dataset)
    clean_anndata_object(anndata_dataset)

    logging.info("Beginning writing to CXG.")
    ctx = tiledb.Ctx(
        {
            "sm.num_reader_threads": 32,
            "sm.num_writer_threads": 32,
            "sm.consolidation.buffer_size": 1 * 1024 * 1024 * 1024,
        }
    )

    tiledb.group_create(cxg_container, ctx=ctx)
    logging.info(f"\t...group created, with name {cxg_container}")

    # Save the dataset metadata
    save_metadata(cxg_container, cxg_group_metadata)
    logging.info("\t...dataset metadata saved")

    # Save the obs/cell dataframe
    save_dataframe(cxg_container, "obs", anndata_dataset.obs, obs_names, ctx=ctx)
    logging.info("\t...obs dataframe created")

    # Save the var/gene dataframe
    save_dataframe(cxg_container, "var", anndata_dataset.var, var_names, ctx=ctx)
    logging.info("\t...var dataframe created")

    # Save the dataset embeddings
    embedding_container = f"{cxg_container}/emb"
    tiledb.group_create(embedding_container, ctx=ctx)
    save_embeddings(embedding_container, anndata_dataset, ctx)
    logging.info("\t...embeddings created")

    # Save the main X matrix of the dataset
    save_x_matrix(cxg_container, anndata_dataset.X, ctx, sparse_threshold)
    logging.info("\t...X created")

    logging.info(f"\t...completed creating CXG container with name {cxg_container}")


def save_metadata(cxg_container, metadata_dict):
    """
    Saves all dataset-wide metadata, including the CXG version and dataset metadata like the title and about link.

    Longer term, tiledb will have support for metadata on groups. Until such feature exists, create an empty array
    and annotate that array.

    https://github.com/TileDB-Inc/TileDB-Py/issues/254
    """

    array_name = f"{cxg_container}/cxg_group_metadata"
    with tiledb.DenseArray(array_name, mode="w") as metadata_array:
        for key, value in metadata_dict.items():
            metadata_array.meta[key] = value


def save_dataframe(container, dataframe_name, dataframe, index_column_name, ctx):
    array_name = f"{container}/{dataframe_name}"

    index_column_name = transform_dataframe_index_into_column(dataframe, dataframe_name, index_column_name)
    create_dataframe(array_name, dataframe)

    with tiledb.DenseArray(array_name, mode="w", ctx=ctx) as array:
        value = {}
        schema_hints = {}
        for key, value in dataframe.items():
            dtype, hints = cxg_type(value)
            value[key] = value.to_numpy(dtype=dtype)
            if hints:
                schema_hints.update({key: hints})

        schema_hints.update({"index": index_column_name})
        array[:] = value
        array.meta["cxg_schema"] = json.dumps(schema_hints)

    tiledb.consolidate(array_name, ctx=ctx)


def save_embeddings(container, anndata_dataset, ctx):
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

    for embedding_name, embedding_values in anndata_dataset.obsm.items():
        if is_valid_embedding(anndata_dataset, embedding_name, embedding_values):
            embedding_name = f"{container}/{embedding_name[2:]}"
            create_embedding(embedding_name, embedding_values)

            with tiledb.DenseArray(embedding_name, mode="w", ctx=ctx) as array:
                array[:] = embedding_values
            tiledb.consolidate(embedding_name, ctx=ctx)
            logging.info(f"\t\t...{embedding_name} embedding created")


def save_x_matrix(container, x_matrix_data, ctx, sparse_threshold):
    """
    Convert and save the main X matrix of the AnnData object.
    """

    x_matrix_name = f"{container}/X"

    shape = x_matrix_data.shape

    # Determine whether the X matrix is sparse or not and if a column shift can be used to make the matrix sparse (or
    # not).
    col_shift = None
    if sparse_threshold == 100.0:
        is_sparse = True
    elif sparse_threshold == 0.0:
        is_sparse = False
    else:
        is_sparse = is_matrix_sparse(x_matrix_data, sparse_threshold)
        if not is_sparse:
            col_shift = is_matrix_sparse_with_column_shift_encoding(x_matrix_data, sparse_threshold)
            is_sparse = col_shift is not None

    create_x(x_matrix_name, shape, is_sparse)
    stride = min(int(np.power(10, np.around(np.log10(1e9 / shape[1])))), 10_000)
    if is_sparse:
        if col_shift is not None:
            logging.info("\t\t...output X as sparse matrix with column shift encoding")
            x_col_shift_name = f"{container}/X_col_shift"
            filters = tiledb.FilterList([tiledb.ZstdFilter()])
            attrs = [tiledb.Attr(dtype=np.float32, filters=filters)]
            domain = tiledb.Domain(tiledb.Dim(domain=(0, shape[1] - 1), tile=min(shape[1], 5000), dtype=np.uint32))
            schema = tiledb.ArraySchema(domain=domain, attrs=attrs)
            tiledb.DenseArray.create(x_col_shift_name, schema)
            with tiledb.DenseArray(x_col_shift_name, mode="w", ctx=ctx) as x_col_shift:
                x_col_shift[:] = col_shift
            tiledb.consolidate(x_col_shift_name, ctx=ctx)
        else:
            logging.info("\t\t...output X as sparse matrix")

        with tiledb.SparseArray(x_matrix_name, mode="w", ctx=ctx) as x_array:
            for start_row_index in range(0, shape[0], stride):
                end_row_index = min(start_row_index + stride, shape[0])
                x_matrix_subset = x_matrix_data[start_row_index:end_row_index, :]
                if not isinstance(x_matrix_subset, np.ndarray):
                    x_matrix_subset = x_matrix_subset.toarray()
                if col_shift is not None:
                    x_matrix_subset = x_matrix_subset - col_shift
                indices = np.nonzero(x_matrix_subset)
                trow = indices[0] + start_row_index
                x_array[trow, indices[1]] = x_matrix_subset[indices[0], indices[1]]

    else:
        logging.info("\t\t...output X as dense matrix")
        with tiledb.DenseArray(x_matrix_name, mode="w", ctx=ctx) as x_array:
            for start_row_index in range(0, shape[0], stride):
                end_row_index = min(start_row_index + stride, shape[0])
                x_matrix_subset = x_matrix_data[start_row_index:end_row_index, :]
                if not isinstance(x_matrix_subset, np.ndarray):
                    x_matrix_subset = x_matrix_subset.toarray()
                x_array[start_row_index:end_row_index, :] = x_matrix_subset

    tiledb.consolidate(x_matrix_name, ctx=ctx)
    if hasattr(tiledb, "vacuum"):
        tiledb.vacuum(x_matrix_name)


def create_x(x_matrix_name, shape, is_sparse):
    """
    The X matrix is accessed in both row and column oriented patterns, depending on the particular operation.
    Because of the data type, default compression works best. The tile size, (50, 100) for dense, and (512,
    2048) for sparse, and global layout (row/col) was chosen empirically, by benchmarking the current cellxgene backend.
    """

    filters = tiledb.FilterList([tiledb.ZstdFilter()])
    attrs = [tiledb.Attr(dtype=np.float32, filters=filters)]
    if is_sparse:
        domain = tiledb.Domain(
            tiledb.Dim(name="obs", domain=(0, shape[0] - 1), tile=min(shape[0], 512), dtype=np.uint32),
            tiledb.Dim(name="var", domain=(0, shape[1] - 1), tile=min(shape[1], 2048), dtype=np.uint32),
        )
    else:
        domain = tiledb.Domain(
            tiledb.Dim(name="obs", domain=(0, shape[0] - 1), tile=min(shape[0], 50), dtype=np.uint32),
            tiledb.Dim(name="var", domain=(0, shape[1] - 1), tile=min(shape[1], 100), dtype=np.uint32),
        )
    schema = tiledb.ArraySchema(
        domain=domain, sparse=is_sparse, attrs=attrs, cell_order="row-major", tile_order="col-major"
    )
    if is_sparse:
        tiledb.SparseArray.create(x_matrix_name, schema)
    else:
        tiledb.DenseArray.create(x_matrix_name, schema)


def is_matrix_sparse(x_matrix_data, sparse_threshold):
    """
    Returns whether `x_matrix_data` is sparse or not (i.e. dense). This is determined by figuring out whether the X
    matrix has a sparsity percentage below the sparse_threshold, returning the number of non-zeros encountered and
    number of elements evaluated.  This function may return before evaluating the whole X matrix if it can be
    determined that X is not sparse enough.
    """

    total_number_of_rows = x_matrix_data.shape[0]
    total_number_of_columns = x_matrix_data.shape[1]
    total_number_of_matrix_elements = total_number_of_rows * total_number_of_columns

    # For efficiency, we count the number of non-zero elements in chunks of the matrix at a time until we hit the
    # maximum number of non zero values allowed before the matrix is deemed "dense." This allows the function the
    # quit early for large dense matrices.
    row_stride = min(int(np.power(10, np.around(np.log10(1e9 / total_number_of_columns)))), 10_000)

    maximum_number_of_non_zero_elements_in_matrix = int(
        total_number_of_rows * total_number_of_columns * sparse_threshold / 100
    )
    number_of_non_zero_elements = 0

    for start_row_index in range(0, total_number_of_rows, row_stride):
        end_row_index = min(start_row_index + row_stride, total_number_of_rows)

        matrix_subset = x_matrix_data[start_row_index:end_row_index, :]
        if not isinstance(matrix_subset, np.ndarray):
            matrix_subset = matrix_subset.toarray()

        number_of_non_zero_elements += np.count_nonzero(matrix_subset)
        if number_of_non_zero_elements > maximum_number_of_non_zero_elements_in_matrix:
            if end_row_index != total_number_of_rows:
                logging.info(
                    "\t\t\t... X matrix is not sparse. Percentage of non-zero elements (estimate): %6.2f"
                    % (100 * number_of_non_zero_elements)
                    / (end_row_index * total_number_of_columns)
                )
            else:
                logging.info(
                    "\t\t\t... X matrix is not sparse. Percentage of non-zero elements (exact): %6.2f"
                    % (100 * number_of_non_zero_elements)
                    / total_number_of_matrix_elements
                )
            return False, number_of_non_zero_elements, end_row_index * total_number_of_columns

    is_sparse = (100.0 * number_of_non_zero_elements / total_number_of_matrix_elements) < sparse_threshold
    return is_sparse


def is_matrix_sparse_with_column_shift_encoding(x_matrix_data, sparse_threshold):
    """Returns a column shift if there is a column shift that allows the given X matrix to be considered as
    sparse. Column shift encoding works by taking the most common value in each column, then subtracting that value from
    each element of the column.  If each column mostly contains its most common value, then the resulting matrix can
    be very sparse.

    This function determines if column shift encoding can be used to transform the X matrix into a sparse matrix with
    a sparsity below the sparse_threshold. If so, returns the array that stores this encoding. This function also
    returns the number of non-zeros encountered and number of elements evaluated.  This function may return before
    evaluating the whole X matrix if it can be determined that X cannot benefit from column shift encoding.
    """

    total_number_of_rows = x_matrix_data.shape[0]
    total_number_of_columns = x_matrix_data.shape[1]
    total_number_of_matrix_elements = total_number_of_rows * total_number_of_columns

    stride = max(1, 128_000_000 // total_number_of_rows)
    column_shift = np.zeros(total_number_of_columns)

    maximum_number_of_non_zero_elements_in_matrix = int(
        total_number_of_rows * total_number_of_columns * sparse_threshold / 100
    )
    number_of_non_zero_elements = 0

    for start_column_index in range(0, total_number_of_columns, stride):
        end_column_index = min(start_column_index + stride, total_number_of_columns)

        matrix_subset = x_matrix_data[:, start_column_index:end_column_index]
        if not isinstance(matrix_subset, np.ndarray):
            matrix_subset = matrix_subset.toarray()

        matrix_subset_mode = mode(matrix_subset)

        column_shift[start_column_index:end_column_index] = matrix_subset_mode.mode
        number_of_non_zero_elements += total_number_of_rows * (end_column_index - start_column_index) - np.sum(
            matrix_subset_mode.count
        )

        if number_of_non_zero_elements > maximum_number_of_non_zero_elements_in_matrix:
            if end_column_index != total_number_of_columns:
                logging.info(
                    "\t\t\t... X matrix is not sparse even with column shift. Percentage of non-zero elements (estimate): %6.2f"
                    % (100 * number_of_non_zero_elements)
                    / (end_column_index * total_number_of_rows)
                )
            else:
                logging.info(
                    "\t\t\t... X matrix is not sparse even with column shift. Percentage of non-zero elements (exact): %6.2f"
                    % (100 * number_of_non_zero_elements)
                    / total_number_of_matrix_elements
                )
            return None

    is_sparse = (100.0 * number_of_non_zero_elements / total_number_of_matrix_elements) < sparse_threshold
    return column_shift if is_sparse else None


def create_dataframe(array_name, dataframe):
    """
    Current access patterns are oriented toward reading very large slices of the dataframe, one attribute at a time.
    Attribute data also tends to be (often) repetitive (bools, categories, strings). Given this, we use a large tile
    size (1000) and very aggressive compression levels.
    """

    tiledb_filter = tiledb.FilterList(
        [
            # Attempt aggressive compression as many of these dataframes are very repetitive strings, bools and other
            # non-float data.
            tiledb.ZstdFilter(level=22),
        ]
    )
    attrs = [
        tiledb.Attr(name=column, dtype=cxg_type(dataframe[column])[0], filters=tiledb_filter) for column in dataframe
    ]
    domain = tiledb.Domain(
        tiledb.Dim(domain=(0, dataframe.shape[0] - 1), tile=min(dataframe.shape[0], 1000), dtype=np.uint32)
    )
    schema = tiledb.ArraySchema(
        domain=domain, sparse=False, attrs=attrs, cell_order="row-major", tile_order="row-major"
    )
    tiledb.DenseArray.create(array_name, schema)


def create_embedding(embedding_name, embedding):
    """
    Embeddings are typically accessed with very large slices (or all of the embedding), and do not benefit from
    overly aggressive compression due to their format.  Given this, we use a large tile size (1000) but only default
    compression level.
    """

    filters = tiledb.FilterList([tiledb.ZstdFilter()])
    attrs = [tiledb.Attr(dtype=embedding.dtype, filters=filters)]
    dimensions = [
        tiledb.Dim(
            domain=(0, embedding.shape[dimension] - 1), tile=min(embedding.shape[dimension], 1000), dtype=np.uint32
        )
        for dimension in range(embedding.ndim)
    ]
    domain = tiledb.Domain(*dimensions)
    schema = tiledb.ArraySchema(
        domain=domain, sparse=False, attrs=attrs, capacity=1_000_000, cell_order="row-major", tile_order="row-major"
    )
    tiledb.DenseArray.create(embedding_name, schema)


def dtype_to_schema(dtype):
    if dtype == np.float32:
        return (np.float32, {})
    if dtype == np.int32:
        return (np.int32, {})
    if dtype == np.bool_:
        return (np.uint8, {"type": "boolean"})
    if dtype == np.str:
        return (np.unicode, {"type": "string"})
    if dtype == "category":
        type = cxg_type(dtype.categories)
        return (type, {"type": "categorical", "categories": dtype.categories.tolist()})

    raise TypeError(f"Annotations of type {dtype} are unsupported.")


def cxg_type(array):
    try:
        return dtype_to_schema(array.dtype)
    except TypeError:
        dtype = array.dtype
        data_kind = dtype.kind
        if array.dtype.kind == "f":
            # Castable to float32
            return np.float32
        if array.dtype.kind in ["i", "u"]:
            # Castable to int32
            if np.can_cast(array.dtype, np.int32):
                return np.int32
            ii32 = np.iinfo(np.int32)
            if array.min() >= ii32.min and array.max() <= ii32.max:
                return np.int32
        if data_kind == "O" and dtype == "object":
            return np.unicode

        raise TypeError(f"Annotations of type {dtype} are unsupported.")


def create_cxg_group_metadata(anndata_dataset, title, about, corpora_properties, extract_colors):
    """
    Return a dictionary containing metadata about CXG dataset. This include data about the version as well as Corpora
    schema properties if they exist, among other pieces of metadata.
    """

    cxg_group_metadata = {"cxg_version": CXG_VERSION, "cxg_properties": json.dumps({"title": title, "about": about})}
    if corpora_properties is not None:
        cxg_group_metadata["corpora"] = json.dumps(corpora_properties)

    if extract_colors:
        try:
            cxg_group_metadata["cxg_category_colors"] = json.dumps(
                convert_anndata_category_colors_to_cxg_category_colors(anndata_dataset)
            )
        except ColorFormatException:
            logging.warning(
                "Failed to extract colors from H5AD file! Fix the H5AD file or rerun with --disable-custom-colors. See help for more details."
            )

    return cxg_group_metadata


def transform_dataframe_index_into_column(dataframe, dataframe_name, index_column_name):
    """
    We rely in the existence of a unique, human-readable index for any dataframe (eg, var is typically gene name,
    obs the cell name). The user can specify these via the --obs-names and --var-names config. If they are not
    specified, use the existing index to create them, giving the resulting column a unique name (e.g. "name").

    In both cases, enforce that the result is unique, and communicate the index column name via the 'index' field in
    the schema hints.
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
        raise KeyError(f"Annotation {index_column_name}, specified in --{dataframe_name}-name, does not exist.")

    return index_column_name


def extract_title_and_about_dataset(title, about, corpora_properties, input_filename):
    """
    The title and about properties of the dataset are set by the following order: if they are explicitly define,
    then use the explicit value. If the dataset is a Corpora-schema based schema, then extract the title and about
    from the corpora_properties. Otherwise, use the input filename (only for title, about will be blank).
    """

    if corpora_properties:
        corpora_project_links = corpora_properties.get("project_links", [])
        corpora_about_link = next(
            (link for link in corpora_project_links if (link.get("link_type", None) == "SUMMARY")), {}
        )
    else:
        corpora_about_link = {}

    filename = path.splitext(path.basename(input_filename))[0]

    title = title if title else corpora_about_link.get("link_name", filename)
    about = about if about else corpora_about_link.get("link_url")

    return title, about


def extract_corpora_properties(anndata_dataset, title, about):
    """
    Extract properties of the anndata object based on the keywords specified in the Corpora schema and return as a
    dictionary.
    """

    corpora_properties = corpora_get_props_from_anndata(anndata_dataset)

    if not corpora_properties:
        # If the return value is None, this means that we were not able to figure out what version of the Corpora
        # schema the object is using and therefore cannot extract any properties.
        raise ValueError("Unknown source file schema version is unsupported.")

    logging.info("Input file appears to be encoded using Corpora schema standards.")
    if title is not None or about is not None:
        logging.warning("Explicit specification of --title or --about will override Corpora schema fields.")

    return corpora_properties


def get_output_directory(input_filename, output_directory, should_overwrite):
    """
    Get the name of the CXG output directory to be created/populated during the dataset conversion.
    """

    if not path.isfile(output_directory) or (path.isfile(output_directory) and should_overwrite):
        if output_directory.endswith(".cxg"):
            return output_directory
        return output_directory + ".cxg"
    if path.isfile(output_directory) and not should_overwrite:
        raise click.BadParameter(
            f"Output directory {output_directory} already exists. If you'd like to overwrite, then run the command "
            f"with the --overwrite flag."
        )

    return path.splitext(input_filename)[1] + ".cxg"


def validate_input_file_type(input_filename):
    """
    Validate that the input file is of a type that we can handle. Currently the only valid file type is `.h5ad`.
    """

    if not input_filename.endswith(".h5ad"):
        raise click.BadParameter(f"Cannot process input file {input_filename}. File must be an H5AD.")


def validate_anndata_object(anndata_dataset):
    if not anndata_dataset.var.index.is_unique:
        raise ValueError("Variable index in AnnData object is not unique - unable to convert.")
    if not anndata_dataset.obs.index.is_unique:
        raise ValueError("Observation index in AnnData object is not unique - unable to convert.")


def sanitize_keys(keys):
    """
    Returns a dictionary mapping of the old keys in the list of `keys` to its new, clean name that is both safe and
    unique.

    We need names to be safe to use as attribute names in tiledb.  See TileDB-Inc/TileDB#1575 and
    TileDB-Inc/TileDB-Py#294.
    TODO: This can be entirely removed once they add proper escaping.
    """

    # Mask out [~/.] and anything outside the ASCII range.
    mask = re.compile(r"[^ -\.0-\[\]-\}]")
    clean_keys_list = [mask.sub("_", key) for key in keys]

    # Dedupe the clean keys list
    deduped_clean_keys_list = []
    for index, clean_key in enumerate(clean_keys_list):
        total_occurrences_of_clean_key = clean_keys_list.count(clean_key)
        total_occurrences_up_until_current_index = clean_keys_list[:index].count(clean_key)
        deduped_clean_keys_list.append(
            clean_key + "_" + str(total_occurrences_up_until_current_index + 1)
            if total_occurrences_of_clean_key > 1
            else clean_key
        )

    return dict(zip(keys, deduped_clean_keys_list))


def sanitize_mapping(mapping):
    """
    Replace keys in the given mapping dictionary with clean and deduped versions.
    """

    clean_keys = sanitize_keys(mapping.keys())
    for old_key, new_key in clean_keys.items():
        if old_key != new_key:
            mapping[new_key] = mapping[old_key]
            del mapping[old_key]


def clean_anndata_object(anndata_dataset):
    """
    Sanitize and dedupe the keys in the obs, var, and obsm attributes of the given AnnData object.
    """

    logging.info("Sanitizing and deduping AnnData object keys.")
    anndata_dataset.obs.rename(columns=sanitize_keys(anndata_dataset.obs.keys().tolist()), inplace=True)
    anndata_dataset.var.rename(columns=sanitize_keys(anndata_dataset.var.keys().tolist()), inplace=True)
    sanitize_mapping(anndata_dataset.obsm)
