import functools
import json
import os.path
import shlex
import shutil
import subprocess
import sys
from subprocess import STDOUT, PIPE
from tempfile import NamedTemporaryFile

import click
import pandas as pd
from click import BadParameter

from server.annotate.annotation_types import AnnotationType
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
    "or S3 path (s3://)",
)
@click.option(
    "-l",
    "--counts-layer",
    help="If specified, raw counts will be read from the AnnData layer of the specified name. If unspecified, "
    "raw counts will be read from `X` matrix, unless 'raw.X' exists, in which case that will be used.",
)
@click.option(
    "-g",
    "--gene-column-name",
    help="The name of the `var` column that contains gene identifiers. The values in this column will be used to match "
    "genes between the query and reference datasets. If not specified, the gene identifiers are expected to exist "
    "in `var.index`.",
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
    "3) `uns` metadata for the prediction operation",
)
@click.option(
    "-n",
    "--run-name",
    type=str,
    help="An optional run name that will be used as a suffix to form the names of new `obs` annotation columns that "
    "will store the predicted annotation values and confidence scores. This can be used to allow multiple "
    "annotation predictions to be run on a single AnnData object.",
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
@click.option(
    "--model-cache-dir",
    default=".models_cache",
    help="Local directory used to store model files that are retrieved from a remote location. Model files will "
    "be read from this directory first, if they exist, to avoid repeating large downloads.",
)
@click.option("--use-model-cache/--no-use-model-cache", default=True)
@click.option(
    "--use-gpu/--no-use-gpu",
    default=True,
    help="Whether to use a GPU for annotation operations (highly recommended, if available).",
)
# TODO: This is a cell type model-specific arg, so not ideal to specify here as a hardcoded option
@click.option(
    "--classifier",
    default="fine",
    help="{fine, coarse, 1, 2, 3, ...}",  # TODO
)
# TODO: This is a cell type model-specific arg, so not ideal to specify here as a hardcoded option
@click.option(
    "--organism",
    type=click.Choice(["Homo sapiens", "Mus musculus"], case_sensitive=True),
    default="Homo sapiens",
    help="",  # TODO
)
@click.help_option("--help", "-h", help="Show this message and exit.")
def annotate(**cli_args):
    _validate_options(cli_args)

    print(f"Reading query dataset {cli_args['input_h5ad_file']}...")

    annotation_prefix = "_".join(
        filter(None, [cli_args.get("annotation_prefix"), cli_args.get("annotation_type"), cli_args.get("run_name")])
    )

    output_h5ad_file = cli_args["input_h5ad_file"] if cli_args["update_h5ad_file"] else cli_args["output_h5ad_file"]

    model_url = cli_args.get("model_url")
    local_model_path = _retrieve_model(cli_args.get("model_cache_dir"), model_url, cli_args.get("use_model_cache"))

    print(f"Annotating {cli_args.get('input_h5ad_file')} with {cli_args.get('annotation_type')}...")

    if cli_args["annotation_type"] == AnnotationType.CELL_TYPE.value:
        predict_args = dict(
            query_dataset_h5ad_path=cli_args.get("input_h5ad_file"),
            output_h5ad_path=output_h5ad_file,
            annotation_prefix=annotation_prefix,
            counts_layer=cli_args.get("counts_layer"),
            gene_column_name=cli_args.get("gene_column_name"),
            classifier=cli_args.get("classifier"),
            organism=cli_args.get("organism"),
            use_gpu=cli_args.get("use_gpu"),
        )
        # Drop args that have values of `None` as these will cause problems when passing into MLflow predict, since it
        # ultimately gets converted into 1-row Pandas DataFrame (None is interpreted as a float type column!)
        predict_args = dict([(k, v) for k, v in predict_args.items() if v is not None])

        # Invoke prediction using MLflow cli, as a separate process.
        # This fully prepares the Python environment that is needed for executing the model.
        # The Python environment will be reused after it is setup once.
        with NamedTemporaryFile(buffering=0) as predict_args_file:
            # write the mlflow predict arguments to a csv file, which will be passed to mlflow cmd
            pd.DataFrame([json.dumps(predict_args)]).to_csv(predict_args_file, index=None)
            predict_args_file.seek(0)

            # run mlflow prediction in subprocess
            predict_cmd = f"mlflow models predict --env-manager virtualenv -m {local_model_path} " \
                          f"-t csv -i {predict_args_file.name}"
            p = subprocess.Popen(
                args=shlex.split(predict_cmd), stdin=predict_args_file, text=True, bufsize=0, stdout=PIPE, stderr=STDOUT
            )

            # display mlflow process output as it runs
            for line in p.stdout:
                print(line.rstrip())

            p.wait()
            if p.returncode == 0:
                print(f"Wrote annotations to {cli_args.get('output_h5ad_file')}")
            else:
                raise f"Annotation failed! exit code={p.returncode}"
    else:
        raise BadParameter(f"unknown annotation type {cli_args['annotation_type']}")


def _retrieve_model(model_cache_dir, model_url, use_cache=True):
    local_cache_model_path = os.path.join(model_cache_dir, os.path.splitext(os.path.basename(model_url))[0])
    if not os.path.exists(local_cache_model_path) or not use_cache:
        print(f"Retrieving model from {model_url}")
        # download from remote source
        with DataLocator(model_url).local_handle() as model_archive_local_path:
            # unpack archive to local cache dir
            shutil.unpack_archive(model_archive_local_path, local_cache_model_path)
    else:
        print(f"Using cached model at {local_cache_model_path}")

    return local_cache_model_path


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
