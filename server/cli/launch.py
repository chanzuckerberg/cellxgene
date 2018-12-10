import logging
from os import devnull
from os.path import splitext, basename
import sys
import warnings
import webbrowser

import click

from server.app.util.errors import ScanpyFileError
from server.app.util.utils import custom_format_warning


@click.command()
@click.argument("data", metavar="<data file>", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--layout", "-l", type=click.Choice(["umap", "tsne"]), default="umap", show_default=True,
              help="Method for layout.")
@click.option("--diffexp", "-d", type=click.Choice(["ttest"]), default="ttest", show_default=True,
              help="Method for differential expression.")
@click.option("--title", "-t", help="Title to display (if omitted will use file name).", metavar="")
@click.option("--verbose", "-v", is_flag=True, default=False, show_default=True,
              help="Provide verbose output, including warnings and all server requests.")
@click.option("--debug", "-d", is_flag=True, default=False, show_default=True,
              help="Run in debug mode.")
@click.option("--open", "-o", "open_browser", is_flag=True, default=False, show_default=True,
              help="Open the web browser after launch.")
@click.option("--port", "-p", help="Port to run server on.", metavar="", default=5005, show_default=True)
@click.option("--obs-names", default=None, metavar="", help="Name of annotation field to use for observations.")
@click.option("--var-names", default=None, metavar="", help="Name of annotation to use for variables.")
@click.option("--host", default="127.0.0.1", help="Host IP address")
@click.option("--max-category-items", default=100, metavar="", show_default=True,
              help="Limits the number of categorical annotation items displayed.")
@click.option("--diffexp-lfc-cutoff", default=0.01, show_default=True,
              help="Relative expression cutoff used when selecting top N differentially expressed genes")
@click.option("--nan-to-num", is_flag=True, default=False, show_default=True,
              help="Replace all floating point NaN with zero, and infinities with finite numbers")
def launch(data, layout, diffexp, title, verbose, debug, obs_names, var_names,
           open_browser, port, host, max_category_items, diffexp_lfc_cutoff,
           nan_to_num):
    """Launch the cellxgene data viewer.
    This web app lets you explore single-cell expression data.
    Data must be in a format that cellxgene expects, read the
    "getting started" guide.

    Examples:

    > cellxgene launch example_dataset/pbmc3k.h5ad --title pbmc3k

    > cellxgene launch <your data file> --title <your title>"""

    # Startup message
    click.echo("[cellxgene] Starting the CLI...")

    # Argument checking
    name, extension = splitext(data)
    if extension != ".h5ad":
        raise click.FileError(basename(data), hint="file type must be .h5ad")

    if debug:
        verbose = True
        open_browser = False
    else:
        warnings.formatwarning = custom_format_warning

    if not verbose:
        sys.tracebacklimit = 0

    if not title:
        file_parts = splitext(basename(data))
        title = file_parts[0]

    # Setup app
    cellxgene_url = f"http://{host}:{port}"
    api_base = f"{cellxgene_url}/api/"

    # Import Flask app
    from server.app.app import app

    app.config.update(
        DATASET_TITLE=title,
        CXG_API_BASE=api_base
    )

    if not verbose:
        log = logging.getLogger("werkzeug")
        log.setLevel(logging.ERROR)

    click.echo(f"[cellxgene] Loading data from {basename(data)}, this may take awhile...")

    # Fix for anaconda python. matplotlib typically expects python to be installed as a framework TKAgg is usually
    # available and fixes this issue. See https://matplotlib.org/faq/virtualenv_faq.html
    import matplotlib as mpl
    mpl.use('TkAgg')
    from server.app.scanpy_engine.scanpy_engine import ScanpyEngine

    args = {
        "layout": layout,
        "diffexp": diffexp,
        "max_category_items": max_category_items,
        "diffexp_lfc_cutoff": diffexp_lfc_cutoff,
        "obs_names": obs_names,
        "var_names": var_names,
        "nan_to_num": nan_to_num
    }

    try:
        app.data = ScanpyEngine(data, args)
    except ScanpyFileError as e:
        raise click.ClickException(f"{e}")

    if open_browser:
        click.echo(f"[cellxgene] Launching! Opening your browser to {cellxgene_url} now.")
        webbrowser.open(cellxgene_url)
    else:
        click.echo(f"[cellxgene] Launching! Please go to {cellxgene_url} in your browser.")

    click.echo("[cellxgene] Type CTRL-C at any time to exit.")

    if not verbose:
        f = open(devnull, 'w')
        sys.stdout = f

    app.run(host=host, debug=debug, port=port, threaded=True)
