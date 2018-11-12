from os.path import expanduser, isdir, isfile, sep, splitext

import click
from numpy import ndarray, unique
from scipy.sparse.csc import csc_matrix


@click.command()
@click.argument("data", nargs=1, metavar="<dataset: file or path to data>", required=True)
@click.option("--layout", "-l", default=["umap", "tsne"], multiple=True, type=click.Choice(["umap", "tsne"]),
              help="Layout algorithm", show_default=True)
@click.option("--recipe", "-r", default="none", type=click.Choice(["none", "seurat", "zheng17"]),
              help="Preprocessing to run.", show_default=True)
@click.option("--output", "-o", default="", help="Save a new file to filename.", metavar="<filename>")
@click.option("--plotting", "-p", default=False, is_flag=True, help="Whether to generate plots.", show_default=True)
@click.option("--sparse", default=False, is_flag=True, help="Whether to force sparsity.", show_default=True)
@click.option("--overwrite", default=False, is_flag=True, help="Allow file overwriting.", show_default=True)
@click.option("--set-obs-names", default="", help="Named field to set as index for obs.", metavar="<name>")
@click.option("--set-var-names", default="", help="Named field to set as index for var.", metavar="<name>")
@click.option("--make-obs-names-unique", default=True, is_flag=True,
              help="Ensure obs index is unique.", show_default=True)
@click.option("--make-var-names-unique", default=True, is_flag=True,
              help="Ensure var index is unique.", show_default=True)
def prepare(data, layout, recipe, output, plotting, sparse, overwrite,
            set_obs_names, set_var_names, make_obs_names_unique, make_var_names_unique):
    """Preprocesses data for use with cellxgene.

    This tool runs a series of scanpy routines for preparing a dataset
    for use with cellxgene. It loads data from different formats
    (h5ad, loom, or a 10x directory), runs dimensionality reduction,
    computes nearest neighbors, computes a layout, performs clustering,
    and saves the results. Includes additional options for naming
    annotations, ensuring sparsity, and plotting results."""

    # collect slow imports here to make CLI startup more responsive
    click.echo("[cellxgene] Starting CLI...")
    import matplotlib
    matplotlib.use("Agg")
    import scanpy.api as sc

    # scanpy settings
    sc.settings.verbosity = 0
    sc.settings.autosave = True

    # check args
    if sparse and not recipe == "none":
        raise click.UsageError("Cannot use a recipe when forcing sparsity")

    output = expanduser(output)

    if not output:
        click.echo("Warning: No file will be saved, to save the results of cellxgene prepare include "
                   "--output <filename> to save output to a new file")
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
            adata.obs_names_make_unique()
        if make_var_names_unique:
            adata.var_names_make_unique()
        if not adata._obs.index.is_unique:
            click.echo("Warning: obs index is not unique")
        if not adata._var.index.is_unique:
            click.echo("Warning: var index is not unique")

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
        try:
            sc.tl.louvain(adata)
        except ModuleNotFoundError:
            click.echo("\nWarning: louvain module is not installed, no clusters will be calculated. "
                       "To fix this please install cellxgene with the optional feature louvain enabled: "
                       "`pip install cellxgene[louvain]`")

    def run_layout(adata):
        if len(unique(adata.obs["louvain"].values)) < 10:
            palette = "tab10"
        else:
            palette = "tab20"

        if "umap" in layout:
            sc.tl.umap(adata)
            if plotting:
                sc.pl.umap(adata, color="louvain", palette=palette, save="_louvain")

        if "tsne" in layout:
            sc.tl.tsne(adata)
            if plotting:
                sc.pl.tsne(adata, color="louvain", palette=palette, save="_louvain")

    def show_step(item):
        names = {
            "make_sparse": "Ensuring sparsity",
            "run_recipe": f"Running preprocessing recipe \"{recipe}\"",
            "run_pca": "Running PCA",
            "run_neighbors": "Calculating neighbors",
            "run_louvain": "Calculating clusters",
            "run_layout": "Computing layout"
        }
        if item is not None:
            return names[item.__name__]

    steps = [make_sparse, run_recipe, run_pca, run_neighbors, run_louvain, run_layout]

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
