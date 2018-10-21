import click

from numpy import unique, ndarray
from scipy.sparse.csc import csc_matrix
from os.path import isfile, isdir, splitext, expanduser

settings = dict(help_option_names=['-h', '--help'])


@click.command()
@click.argument('dataset', nargs=1, metavar='<dataset: file or path to data>', required=True)
@click.option('--layout', default='umap', type=click.Choice(['umap', 'tsne', 'umap+tsne']),
              help='layout algorithm', show_default=True)
@click.option('--recipe', default='none', type=click.Choice(['none', 'seurat', 'zheng17']),
              help='preprocessing to run', show_default=True)
@click.option('--output', default='', help='save a new file to filename')
@click.option('--sparse/--no-sparse', default=False, help='whether to force sparsity', show_default=True)
@click.option('--overwrite/--no-overwrite', default=False, help='allow overwriting an existing file', show_default=True)
@click.option('--plotting/--no-plotting', default=False, help='whether to generate plots', show_default=True)
def cli(dataset, layout, recipe, output, sparse, overwrite, plotting):
    """
    preprocesses data for use with cellxgene
    """

    # collect slow imports here to make CLI startup more responsive
    click.echo('[cellxgene] Starting CLI...')
    import matplotlib
    matplotlib.use('Agg')
    import scanpy.api as sc

    # scanpy settings
    sc.settings.verbosity = 0
    sc.settings.autosave = True

    # check args
    if sparse and not recipe == 'none':
        raise click.UsageError('Cannot use a recipe when forcing sparsity')

    output = expanduser(output)
    if isfile(output) and not overwrite:
        raise click.UsageError('Cannot overwwrite existing file %s, try using the flag --overwrite' % output)

    def load_data(dataset):
        if isfile(dataset):
            name, extension = splitext(dataset)
            if extension == '.h5ad':
                adata = sc.read_h5ad(dataset)
            elif extension == '.loom':
                adata = sc.read_loom(dataset)
            else:
                raise click.FileError(dataset, hint='does not have a valid extension [.h5ad | .loom]')
        elif isdir(dataset):
            if not dataset.endswith('/'):
                raise click.FileError(dataset, hint='a path must end with a /')
            else:
                adata = sc.read_10x_mtx(dataset)
        else:
            raise click.FileError(dataset, hint='not a valid file or path')

        adata.var_names_make_unique()

        return adata

    def make_sparse(adata):
        if (type(adata.X) is ndarray) and sparse:
            adata.X = csc_matrix(adata.X)

    def run_recipe(adata):
        if recipe == 'seurat':
            sc.pp.recipe_seurat(adata)
        elif recipe == 'zheng17':
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
            sc.pp.pca(adata, svd_solver='arpack', zero_center=False)
        else:
            sc.pp.pca(adata, svd_solver='arpack')

    def run_neighbors(adata):
        sc.pp.neighbors(adata)

    def run_louvain(adata):
        sc.tl.louvain(adata)

    def run_layout(adata):
        if len(unique(adata.obs['louvain'].values)) < 10:
            palette = 'tab10'
        else:
            palette = 'tab20'

        if layout == 'umap' or layout == 'umap+tsne':
            sc.tl.umap(adata)
            if plotting:
                sc.pl.umap(adata, color='louvain', palette=palette, save='_louvain')

        if layout == 'tsne' or layout == 'umap+tsne':
            sc.tl.tsne(adata)
            if plotting:
                sc.pl.tsne(adata, color='louvain', palette=palette, save='_louvain')

    def show_step(item):
        names = {
            'make_sparse': 'Ensuring sparsity',
            'run_recipe': 'Running preprocessing recipe "%s"' % recipe,
            'run_pca': 'Running PCA',
            'run_neighbors': 'Calculating neighbors',
            'run_louvain': 'Calculating clusters',
            'run_layout': 'Computing layout'
        }
        if item is not None:
            return names[item.__name__]

    steps = [make_sparse, run_recipe, run_pca, run_neighbors, run_louvain, run_layout]

    click.echo('[cellxgene] Loading data from %s, please wait...' % dataset)
    adata = load_data(dataset)

    click.echo('[cellxgene] Beginning preprocessing...')
    with click.progressbar(steps, label='[cellxgene] Progress', show_eta=False, item_show_func=show_step) as bar:
        for step in bar:
            step(adata)

    # saving
    if not output == '':
        click.echo('[cellxgene] Saving results to %s...' % output)
        adata.write(output)

    click.echo('[cellxgene] ' + click.style('Success!', fg='green'))


if __name__ == '__main__':
    cli()
