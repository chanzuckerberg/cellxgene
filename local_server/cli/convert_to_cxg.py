from os import path

import click

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
    "input-file", nargs=1, type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "-o",
    "--output-directory",
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
@click.option(
    "--obs-names",
    help="Name to a column in the obs dataframe that will be used as the index for the dataframe instead of "
    "the one designated by the dataframe generated-index.",
)
@click.option(
    "--var-names",
    help="Name to a column in the var dataframe that will be used as the index for the dataframe instead of "
    "the one designated by the dataframe generated-index.",
)
@click.option(
    "--disable-custom-colors",
    help="When set, conversion process will not extract scanpy-compatible category colors from the H5AD file.",
    default=False,
    show_default=True,
    is_flag=True,
)
@click.option(
    "--disable-corpora-schema",
    help="When set, conversion process will neither extract nor store Corpora schema information. See "
    "https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema.md for more "
    "information.",
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
    overwrite,
):
    """
    Convert a dataset file into CXG.
    """

    h5ad_data_file = H5ADDataFile(
        input_file, backed, title, about, obs_names, var_names, use_corpora_schema=not disable_corpora_schema
    )

    # Get the directory that will hold all the CXG files
    cxg_output_container = get_output_directory(input_file, output_directory, overwrite)

    h5ad_data_file.to_cxg(
        cxg_output_container, sparse_threshold, convert_anndata_colors_to_cxg_colors=not disable_custom_colors
    )


def get_output_directory(input_filename, output_directory, should_overwrite):
    """
    Get the name of the CXG output directory to be created/populated during the dataset conversion.
    """

    if output_directory and (not path.isdir(output_directory) or (path.isdir(output_directory) and should_overwrite)):
        if output_directory.endswith(".cxg"):
            return output_directory
        return output_directory + ".cxg"
    if output_directory and path.isdir(output_directory) and not should_overwrite:
        raise click.BadParameter(
            f"Output directory {output_directory} already exists. If you'd like to overwrite, then run the command "
            f"with the --overwrite flag."
        )

    return path.splitext(input_filename)[0] + ".cxg"
