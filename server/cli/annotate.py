import functools
import pickle
import sys
from typing import Any

import click
import scanpy as sc
from anndata import AnnData
from click import BadParameter
from numpy import unique

from server.annotate.annotation_types import AnnotationType

from server.annotate import cell_type
from server.common.utils.data_locator import DataLocator

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
    required=True,
    type=str,
    help="The input H5AD file containing the missing annotations.",
)
@click.option(
    "-m",
    "--model-url",
    required=True,
    help="The URL of the model used to prediction annotated labels. May be a local filesystem directory "
         "or S3 path (s3://)"
)
@click.option(
    "-l",
    "--h5ad-layer",
    help="If specified, raw counts will be read from the AnnData layer of the specified name. If unspecified, "
         "raw counts will be read from `X` matrix, unless 'raw.X' exists, in which case that will be used."
)
@click.option(
    "-g",
    "--gene-column-name",
    help="The name of the `var` column that contains gene identifiers. The values in this column will be used to match "
         "genes between the query and reference datasets. If not specified, the gene identifiers are expected to exist "
         "in `var.index`."
)
# TODO: Useful if we want to support discoverability of models
# @click.option(
#     "-r",
#     "--model-repository",
#     help="The base URL of the model repository. Maybe a local filesystem directory or S3 path (s3://)"
# )
# TODO: Useful if we want to support other, future annotation types, beyond "Cell Type". Currently hidden
@click.option(
    "-a",
    "--annotation-type",
    type=click.Choice([t.value for t in AnnotationType]),
    default=AnnotationType.CELL_TYPE.value,
    show_default=True,
    hidden=True,  # Remove if we add support for more annotation types
    help="The type of annotation to perform. This model to be used will be inferred from the annotation type.",
)
@click.option(
    "-c",
    "--annotation-column-prefix",
    type=str,
    default="cxg_predicted",
    show_default=True,
    help="An optional prefix used to form the names of new `obs` annotation columns that will store the predicted "
         "annotation values and confidence scores."
)
@click.option(
    "-n",
    "--run-name",
    type=str,
    help="An optional run name that will be used as a suffix to form the names of new `obs` annotation columns that "
         "will store the predicted annotation values and confidence scores. This can be used to allow multiple "
         "annotation predictions to be run on a single AnnData object."
)
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

    _validate_options(cli_args)

    query_dataset = sc.read_h5ad(cli_args['input_h5ad_file'])
    model = _fetch_model(cli_args['model_url'])
    annotation_column_name = '_'.join(filter(None,
                                             [cli_args.get('annotation_column_prefix'),
                                              cli_args.get('annotation_type'),
                                              cli_args.get('run_name')]))

    if cli_args['annotation_type'] == AnnotationType.CELL_TYPE.value:
        cell_type.annotate(query_dataset, model, annotation_column_name=annotation_column_name,
                           gene_col_name=cli_args['gene_column_name'],
                           min_common_gene_count=1)  # TODO
    else:
        raise BadParameter(f"unknown annotation type {cli_args['annotation_type']}")

    output_h5ad_file = cli_args['input_h5ad_file'] if cli_args['update_h5ad_file'] else cli_args['output_h5ad_file']
    query_dataset.write_h5ad(output_h5ad_file)

    print(f"added annotations to {output_h5ad_file} in obs column '{annotation_column_name}'")


def _validate_options(cli_args):
    # TODO(atolopko): Use cloup library for this logic
    if cli_args["update_h5ad_file"] and cli_args["output_h5ad_file"]:
        click.echo("--update_h5ad_file and --output_h5ad_file are mutually exclusive")
        sys.exit(1)
    if not (cli_args["update_h5ad_file"] or cli_args["output_h5ad_file"]):
        click.echo("--update_h5ad_file or --output_h5ad_file must be specified")
        sys.exit(1)


def _fetch_model(model_url) -> Any:
    model_url_locator = DataLocator(model_url)
    if not model_url_locator.exists():
        raise f"model file '{model_url}' not found"

    with model_url_locator.open() as f:
        return pickle.load(f)


if __name__ == "__main__":
    annotate()
