import functools
import sys

import click

from server.common.utils.utils import sort_options


def annotate_args(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


@sort_options
@click.command(
    short_help="Annotate H5AD file columns. Run `cellxgene annotation --help` for more information.",
    options_metavar="<options>",
)
@click.option(
    "-i",
    "--input-h5ad-file",
    type=str,
    help="The input H5AD file containing the missing annotations.",
)
@click.option(
    "-m",
    "--model-url",
    help="The URL of the model used to prediction annotated labels. May be a local filesystem directory "
         "or S3 path (s3://)"
)
@click.option(
    "-l",
    "--h5ad-layer",
    help="If specified, raw counts will be read from the AnnData layer of the specified name. If unspecified, "
         "raw counts will be read from `X` matrix, unless 'raw.X' exists, in which case that will be used."

)
# TODO: Useful if we want to support multiple tissue types (for cell type annotation models)
# @click.option(
#     "-r",
#     "--model-repository",
#     help="The base URL of the model repository. Maybe a local filesystem directory or S3 path (s3://)"
# )
# TODO: Useful if we want to support other, future annotation types, beyond "Cell Type"
# @click.option(
#     "-a",
#     "--annotation-type",
#     type=click.Choice([t.value for t in AnnotationType]),
#     default=AnnotationType.CELL_TYPE.value,
#     help="The type of annotation to perform. This model to be used will be inferred from the annotation type.",
# )
@click.option(
    "-c",
    "--annotation-column-prefix",
    type=str,
    default="cxg_predicted_cell_type_",
    help="A prefix used to form the names of new `obs` annotation columns that will store the predicted annotation "
         "values and confidence scores."
)
# @click.option(
#     "-s",
#     "--annotation-column-suffix",
#     type=str,
#     help="An optional suffix used to form the names of new `obs` annotation columns that will store the predicted "
#          "annotation values and confidence scores. This can be used to allow multiple annotation predictions to be "
#          "run on a single AnnData object."
# )
@click.option(
    "-u",
    "--update-h5ad-file",
    is_flag=True,
    help="Flag indicating whether to update the input h5ad file with annotation values.  This option is mutually "
    "exclusive with --output-h5ad-file.",
)
@click.option(
    "-o",
    "--output-h5ad-file",
    help="The output H5AD file that will contain the generated annotation values. This option is mutually "
    "exclusive with --update-h5ad-file.",
)
@click.help_option("--help", "-h", help="Show this message and exit.")
def annotate(**cli_args):
    print(cli_args)

    # TODO(atolopko): Use cloup library for this logic
    if cli_args["update_h5ad_file"] and cli_args["output_h5ad_file"]:
        click.echo("--update_h5ad_file and --output_h5ad_file are mutually exclusive")
        sys.exit(1)
    if not (cli_args["update_h5ad_file"] or cli_args["output_h5ad_file"]):
        click.echo("--update_h5ad_file or --output_h5ad_file must be specified")
        sys.exit(1)


if __name__ == "__main__":
    annotate()
