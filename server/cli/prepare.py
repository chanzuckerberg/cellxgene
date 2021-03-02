from os.path import expanduser, isdir, isfile, sep, splitext

import click
import pandas as pd
from numpy import ndarray, unique
from scipy.sparse.csc import csc_matrix

from shared_utils.utils.utils import sort_options


@sort_options
@click.command(
    short_help="Preprocess data for use with cellxgene. " "Run `cellxgene prepare --help` for more information.",
    options_metavar="<options>",
)
@click.argument("data", nargs=1, metavar="<path to data file>", required=True)
@click.option(
    "--embedding",
    "-e",
    default=["umap", "tsne"],
    multiple=True,
    type=click.Choice(["umap", "tsne"]),
    help="Embedding algorithm(s). Repeat option for multiple embeddings.",
    show_default=True,
)
@click.option(
    "--recipe", "-r", default="none", type=click.Choice(["none", "seurat", "zheng17"]), show_default=True,
)
@click.option("--output", "-o", default="", help="Save a new file to filename.", metavar="<filename>")
@click.option("--plotting", "-p", default=False, is_flag=True, help="Generate plots.", show_default=True)
@click.option("--sparse", default=False, is_flag=True, help="Force sparsity.", show_default=True)
@click.option("--overwrite", default=False, is_flag=True, help="Allow file overwriting.", show_default=True)
@click.option("--set-obs-names", default="", help="Named field to set as index for obs.", metavar="<name>")
@click.option("--set-var-names", default="", help="Named field to set as index for var.", metavar="<name>")
@click.option(
    "--skip-qc",
    default=False,
    is_flag=True,
    help="Do not run quality control metrics. By default cellxgene runs them "
    "(saved to adata.obs and adata.var; see scanpy.pp.calculate_qc_metrics for details).",
)
@click.option(
    "--make-obs-names-unique/--no-make-obs-names-unique",
    default=True,
    help="Ensure obs index is unique.",
    show_default=True,
)
@click.option(
    "--make-var-names-unique/--no-make-var-names-unique",
    default=True,
    help="Ensure var index is unique.",
    show_default=True,
)
@click.help_option("--help", "-h", help="Show this message and exit.")
def prepare(
    data,
    embedding,
    recipe,
    output,
    plotting,
    sparse,
    overwrite,
    set_obs_names,
    set_var_names,
    skip_qc,
    make_obs_names_unique,
    make_var_names_unique,
):
    """
    Preprocess data for use with cellxgene.
    This tool runs a series of scanpy routines for preparing a dataset for use
    with cellxgene. It loads data from different formats
    (h5ad, loom, or a 10x directory), runs dimensionality reduction,
    computes nearest neighbors, computes an embedding, performs clustering,
    and saves the results. Includes additional options for naming annotations,
    ensuring sparsity, and plotting results.
    """

    # collect slow imports here to make CLI startup more responsive
    click.echo("[cellxgene] Starting CLI...")
    try:
        import matplotlib

        matplotlib.use("Agg")
        import scanpy as sc
    except ImportError:
        raise click.ClickException(
            "[cellxgene] cellxgene prepare has not been installed. Please run `pip install 'cellxgene[prepare]'` "
            "to install the necessary requirements."
        )

    # scanpy settings
    sc.settings.verbosity = 0
    sc.settings.autosave = True

    # check args
    if sparse and not recipe == "none":
        raise click.UsageError("Cannot use a recipe when forcing sparsity")

    output = expanduser(output)

    if not output:
        click.echo(
            "Warning: No file will be saved, to save the results of cellxgene prepare include "
            "--output <filename> to save output to a new file"
        )
    if isfile(output) and not overwrite:
        raise click.UsageError(f"Cannot overwrite existing file {output}, try using the flag --overwrite")

    def load_data(data):
        if isfile(data):
            name, extension = splitext(data)
            if extension == ".h5ad":
                adata = sc.read_h5ad(data)
            elif extension == ".loom":
                adata = sc.read_loom(data)
            else:
                raise click.FileError(data, hint="does not have a valid extension [.h5ad | .loom]")
        elif isdir(data):
            if not data.endswith(sep):
                data += sep
            adata = sc.read_10x_mtx(data)
        else:
            raise click.FileError(data, hint="not a valid file or path")

        if not set_obs_names == "":
            if set_obs_names not in adata.obs_keys():
                raise click.UsageError(f"obs {set_obs_names} not found, options are: {adata.obs_keys()}")
            adata.obs_names = adata.obs[set_obs_names]
        if not set_var_names == "":
            if set_var_names not in adata.var_keys():
                raise click.UsageError(f"var {set_var_names} not found, options are: {adata.var_keys()}")
            adata.var_names = adata.var[set_var_names]
        if make_obs_names_unique:
            adata.obs.index = make_index_unique(adata.obs.index)
        if make_var_names_unique:
            adata.var.index = make_index_unique(adata.var.index)
        if not adata._obs.index.is_unique:
            click.echo("Warning: obs index is not unique")
        if not adata._var.index.is_unique:
            click.echo("Warning: var index is not unique")
        return adata

    def calculate_qc_metrics(adata):
        if not skip_qc:
            sc.pp.calculate_qc_metrics(adata, inplace=True)
        return adata

    def make_sparse(adata):
        if (type(adata.X) is ndarray) and sparse:
            adata.X = csc_matrix(adata.X)

    def run_recipe(adata):
        if recipe == "seurat":
            sc.pp.recipe_seurat(adata)
        elif recipe == "zheng17":
            sc.pp.recipe_zheng17(adata)
        else:
            sc.pp.filter_cells(adata, min_genes=5)
            sc.pp.filter_genes(adata, min_cells=25)
            if sparse:
                sc.pp.scale(adata, zero_center=False)
            else:
                sc.pp.scale(adata)

    def run_pca(adata):
        if sparse:
            sc.pp.pca(adata, svd_solver="arpack", zero_center=False)
        else:
            sc.pp.pca(adata, svd_solver="arpack")

    def run_neighbors(adata):
        sc.pp.neighbors(adata)

    def run_louvain(adata):
        sc.tl.louvain(adata)

    def run_embedding(adata):
        if len(unique(adata.obs["louvain"].values)) < 10:
            palette = "tab10"
        else:
            palette = "tab20"

        if "umap" in embedding:
            sc.tl.umap(adata)
            if plotting:
                sc.pl.umap(adata, color="louvain", palette=palette, save="_louvain")

        if "tsne" in embedding:
            sc.tl.tsne(adata)
            if plotting:
                sc.pl.tsne(adata, color="louvain", palette=palette, save="_louvain")

    def show_step(item):
        if not skip_qc:
            qc_name = "Calculating QC metrics"
        else:
            qc_name = "Skipping QC"
        names = {
            "calculate_qc_metrics": qc_name,
            "make_sparse": "Ensuring sparsity",
            "run_recipe": f'Running preprocessing recipe "{recipe}"',
            "run_pca": "Running PCA",
            "run_neighbors": "Calculating neighbors",
            "run_louvain": "Calculating clusters",
            "run_embedding": "Computing embedding",
        }
        if item is not None:
            return names[item.__name__]

    steps = [calculate_qc_metrics, make_sparse, run_recipe, run_pca, run_neighbors, run_louvain, run_embedding]

    click.echo(f"[cellxgene] Loading data from {data}, please wait...")
    adata = load_data(data)

    click.echo("[cellxgene] Beginning preprocessing...")
    with click.progressbar(steps, label="[cellxgene] Progress", show_eta=False, item_show_func=show_step) as bar:
        for step in bar:
            step(adata)

    # saving
    if not output == "":
        click.echo(f"[cellxgene] Saving results to {output}...")
        adata.write(output)

    click.echo("[cellxgene] Success!")


# TODO (mweiden): remove this once this issue is resolved https://github.com/theislab/anndata/issues/344
# Note: tentative solution here https://github.com/theislab/anndata/pull/345
def make_index_unique(index: pd.Index, join: str = "-"):
    """
    Makes the index unique by appending a number string to each duplicate index element: '1', '2', etc.

    If a tentative name created by the algorithm already exists in the index, it tries the next integer in the sequence.

    The first occurrence of a non-unique value is ignored.
    Parameters
    ----------
    join
         The connecting string between name and integer.
    Examples
    --------
    >>> from anndata import AnnData
    >>> adata1 = AnnData(np.ones((3, 2)), dict(obs_names=['a', 'b', 'c']))
    >>> adata2 = AnnData(np.zeros((3, 2)), dict(obs_names=['d', 'b', 'b']))
    >>> adata = adata1.concatenate(adata2)
    >>> adata.obs_names
    Index(['a', 'b', 'c', 'd', 'b', 'b'], dtype='object')
    >>> adata.obs_names_make_unique()
    >>> adata.obs_names
    Index(['a', 'b', 'c', 'd', 'b-1', 'b-2'], dtype='object')
    """
    if index.is_unique:
        return index
    from collections import defaultdict

    values = index.values
    values_set = set(values)
    indices_dup = index.duplicated(keep="first")
    values_dup = values[indices_dup]
    counter = defaultdict(lambda: 0)
    for i, v in enumerate(values_dup):
        while True:
            counter[v] += 1
            tentative_new_name = v + join + str(counter[v])
            if tentative_new_name not in values_set:
                values_set.add(tentative_new_name)
                values_dup[i] = tentative_new_name
                break

    values[indices_dup] = values_dup
    index = pd.Index(values)
    return index
