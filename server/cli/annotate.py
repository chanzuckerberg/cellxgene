import functools
import json
import sys

import click
import scanpy as sc
from click import BadParameter

from server.annotate import cell_type
from server.annotate.annotation_types import AnnotationType
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
    "--counts-layer",
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
    "--annotation-prefix",
    type=str,
    default="cxg",
    show_default=True,
    help="An optional prefix used to form the names of: 1) new `obs` annotation columns that will store the predicted "
         "annotation values and confidence scores, 2) `obsm` embeddings (reference and umap embedding), and "
         "3) `uns` metadata for the prediction operation"
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
# This is fixed by scArches: scvi/model/base/_archesmixin.py:203
# @click.option(
#     "--min-common-gene-pct",
#     default=80,
#     type=click.IntRange(0, 100, min_open=True),
#     show_default=True,
#     help="The percentage of reference dataset genes that must also exist in the query dataset. This command will abort "
#          "if there are fewer genes found to be in common. Common genes will be determined by matching the var.index "
#          "column values of the reference and query datasets, unless the --gene-column-name option is specified, in "
#          "which case the reference dataset's var.index will matched to the specified column in the query dataset's "
#          "`obs.var` axis. exclusive with --update-h5ad-file."
# )
@click.option(
    "--model-cache-dir",
    default=".models_cache",
    help="Local directory used to store model files that are retrieved from a remote location. Model files will "
         "be read from this directory first, if they exist, to avoid repeating large downloads."
)
@click.option(
    "--use-gpu/--no-use-gpu",
    default=True,
    help="Whether to use a GPU for annotation operations (highly recommended, if available)."
)
@click.option(
    "--scarches-train-param-overrides",
    default='{}',
    hidden=True
)
@click.help_option("--help", "-h", help="Show this message and exit.")
def annotate(**cli_args):
    _validate_options(cli_args)

    print(f"Reading query dataset {cli_args['input_h5ad_file']}...")
    query_dataset = sc.read_h5ad(cli_args['input_h5ad_file'])

    annotation_prefix = '_'.join(filter(None,
                                        [cli_args.get('annotation_prefix'),
                                         cli_args.get('annotation_type'),
                                         cli_args.get('run_name')]))

    model_url = cli_args.get('model_url')
    model_cache_dir = cli_args.get('model_cache_dir')

    print(f"Annotating {cli_args.get('input_h5ad_file')} with {cli_args.get('annotation_type')}...")

    if cli_args['annotation_type'] == AnnotationType.CELL_TYPE.value:
        cell_type.annotate(query_dataset=query_dataset,
                           model_url=model_url,
                           model_cache_dir=model_cache_dir,
                           annotation_prefix=annotation_prefix,
                           gene_column_name=cli_args.get('gene_column_name'),
                           counts_layer=cli_args.get('counts_layer'),
                           use_gpu=cli_args['use_gpu'],
                           train_param_overrides=json.loads(cli_args['scarches_train_param_overrides'])
                           )
    else:
        raise BadParameter(f"unknown annotation type {cli_args['annotation_type']}")

    output_h5ad_file = cli_args['input_h5ad_file'] if cli_args['update_h5ad_file'] else cli_args['output_h5ad_file']
    query_dataset.write_h5ad(output_h5ad_file)


def _validate_options(cli_args):
    # TODO(atolopko): Use cloup library for this logic
    if cli_args["update_h5ad_file"] and cli_args["output_h5ad_file"]:
        click.echo("--update_h5ad_file and --output_h5ad_file are mutually exclusive")
        sys.exit(1)
    if not (cli_args["update_h5ad_file"] or cli_args["output_h5ad_file"]):
        click.echo("--update_h5ad_file or --output_h5ad_file must be specified")
        sys.exit(1)


if __name__ == "__main__":
    annotate()
