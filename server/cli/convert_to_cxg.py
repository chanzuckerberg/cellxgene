from os import path

import click
import numpy as np

from server.converters.h5ad_data_file import H5ADDataFile


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

    h5ad_data_file = H5ADDataFile(input_file, backed, title, about, obs_names, var_names, not disable_corpora_schema)

    # Get the directory that will hold all the CXG files
    cxg_output_container = get_output_directory(input_file, output_directory, should_overwrite)

    h5ad_data_file.to_cxg(cxg_output_container, not disable_custom_colors, sparse_threshold)


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
