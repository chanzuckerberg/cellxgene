import click
import logging
import webbrowser

from os.path import splitext, basename, isfile


@click.command()
@click.argument('data', metavar='<data file>')
@click.option('--layout', '-l', type=click.Choice(['umap', 'tsne']), default='umap', show_default=True,
              help='Method for layout.')
@click.option('--diffexp', '-d', type=click.Choice(['ttest']), default='ttest', show_default=True,
              help='Method for differential expression.')
@click.option('--title', '-t', help='Title to display (if omitted will use file name).', metavar='')
@click.option('--verbose', '-v', is_flag=True, default=False, show_default=True,
              help='Provide verbose output, including warnings and all server requests.')
@click.option('--debug', '-d', is_flag=True, default=False, show_default=True,
              help='Run in debug mode.')
@click.option('--open', '-o', 'open_browser', is_flag=True, default=False, show_default=True,
              help='Open the web browser after launch.')
@click.option('--port', '-p', help="Port to run server on.", metavar='', default=5005, show_default=True)
@click.option('--listen-all', is_flag=True, default=False, show_default=True,
              help='Bind to all interfaces (this makes the server accessible beyond this computer).')
@click.option('--max-category-items', default=100, metavar='', show_default=True,
              help='Limits the number of categorical annotation items displayed.')
def launch(data, layout, diffexp, title, verbose, debug, open_browser, port, listen_all, max_category_items):
    """
    Launch the cellxgene data viewer.
    This web app lets you explore single-cell expression data.
    Data must be in a format that cellxgene expects, read the
    "getting started" guide.

    Examples:

    > cellxgene launch example_dataset/pbmc3k.h5ad --title pbmc3k

    > cellxgene launch <your data file> --title <your title>
    """
    # Startup message
    click.echo('[cellxgene] Starting the CLI...')

    # Import Flask app
    from ..app.app import app

    # Argument checking
    if not isfile(data):
        raise click.FileError(data, hint='file does not exist')
    else:
        name, extension = splitext(data)
        if not extension == '.h5ad':
            raise click.FileError(data, hint='file type must be .h5ad')

    if debug:
        verbose = True
        open_browser = False

    if not title:
        file_parts = splitext(basename(data))
        title = file_parts[0]

    if listen_all:
        host = '0.0.0.0'
    else:
        host = '127.0.0.1'

    # Setup app
    cellxgene_url = f"http://{host}:{port}"
    api_base = f"{cellxgene_url}/api/"

    app.config.update(
        DATASET_TITLE=title,
        CXG_API_BASE=api_base
        )

    if not verbose:
        log = logging.getLogger('werkzeug')
        log.setLevel(logging.ERROR)

    click.echo('[cellxgene] Loading data from %s, this may take awhile...' % data)

    from ..app.scanpy_engine.scanpy_engine import ScanpyEngine
    app.data = ScanpyEngine(data, layout_method=layout, diffexp_method=diffexp,
                            max_category_items=max_category_items)

    if open_browser:
        click.echo('[cellxgene] Launching! Opening your browser to %s now.' % cellxgene_url)
        webbrowser.open(cellxgene_url)
    else:
        click.echo('[cellxgene] Launching! Please go to %s in your browser.' % cellxgene_url)

    click.echo('[cellxgene] Type CTRL-C at any time to exit.')

    app.run(host=host, debug=debug, port=port)
