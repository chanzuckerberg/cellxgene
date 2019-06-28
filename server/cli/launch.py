import errno
import functools
import logging
from os import devnull
from os.path import splitext, basename, getsize
import sys
import warnings
import webbrowser

import click

from server.app.app import Server
from server.app.util.errors import ScanpyFileError
from server.app.util.utils import custom_format_warning
from server.utils.utils import find_available_port, is_port_available

# anything bigger than this will generate a special message
BIG_FILE_SIZE_THRESHOLD = 100 * 2 ** 20  # 100MB


def common_args(func):
    """
    Decorator to contain CLI args that will be common to both CLI and GUI: title and engine args.
    """
    @click.option("--title", "-t", help="Title to display (if omitted will use file name).")
    @click.option(
        "--layout",
        "-l",
        default=[],
        multiple=True,
        show_default=True,
        help="Layout name, eg, 'umap'."
    )
    @click.option("--obs-names", default=None, metavar="", help="Name of annotation field to use for observations.")
    @click.option("--var-names", default=None, metavar="", help="Name of annotation to use for variables.")
    @click.option(
        "--max-category-items",
        default=1000,
        metavar="",
        show_default=True,
        help="Categories with more distinct values than this will not be displayed.",
    )
    @click.option(
        "--diffexp-lfc-cutoff",
        default=0.01,
        show_default=True,
        help="Relative expression cutoff used when selecting top N differentially expressed genes",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper


def parse_engine_args(layout, obs_names, var_names, max_category_items, diffexp_lfc_cutoff):
    return {
        "layout": layout,
        "max_category_items": max_category_items,
        "diffexp_lfc_cutoff": diffexp_lfc_cutoff,
        "obs_names": obs_names,
        "var_names": var_names,
    }


@click.command()
@click.argument("data", metavar="<data file>", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    show_default=True,
    help="Provide verbose output, including warnings and all server requests.",
)
@click.option("--debug", is_flag=True, default=False, show_default=True, help="Run in debug mode.")
@click.option(
    "--open",
    "-o",
    "open_browser",
    is_flag=True,
    default=False,
    show_default=True,
    help="Open the web browser after launch.",
)
@click.option("--port", "-p", help="Port to run server on, if not specified cellxgene will find an available port.",
              metavar="", show_default=True)
@click.option("--host", default="127.0.0.1", help="Host IP address")
@click.option(
    "--scripts",
    default=[],
    multiple=True,
    help="Additional script files to include in html page",
    show_default=True,
)
@common_args
def launch(
        data,
        verbose,
        debug,
        open_browser,
        port,
        host,
        layout,
        obs_names,
        var_names,
        max_category_items,
        diffexp_lfc_cutoff,
        title,
        scripts
):
    """Launch the cellxgene data viewer.
    This web app lets you explore single-cell expression data.
    Data must be in a format that cellxgene expects, read the
    "getting started" guide.

    Examples:

    > cellxgene launch example_dataset/pbmc3k.h5ad --title pbmc3k

    > cellxgene launch <your data file> --title <your title>"""

    e_args = parse_engine_args(layout, obs_names, var_names, max_category_items, diffexp_lfc_cutoff)
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

    if scripts:
        click.echo(r"""
    / / /\ \ \__ _ _ __ _ __ (_)_ __   __ _
    \ \/  \/ / _` | '__| '_ \| | '_ \ / _` |
     \  /\  / (_| | |  | | | | | | | | (_| |
      \/  \/ \__,_|_|  |_| |_|_|_| |_|\__, |
                                      |___/
    The --scripts flag is intended for developers to include google analytics etc. You could be opening yourself to a
    security risk by including the --scripts flag. Make sure you trust the scripts that you are including.
            """)
        scripts_pretty = ", ".join(scripts)
        click.confirm(f"Are you sure you want to inject these scripts: {scripts_pretty}?", abort=True)

    if not title:
        file_parts = splitext(basename(data))
        title = file_parts[0]

    if port:
        if debug:
            raise click.ClickException("--port and --debug may not be used together (try --verbose for error logging).")
        if not is_port_available(host, int(port)):
            raise click.ClickException(
                f"The port selected {port} is in use, please specify an open port using the --port flag."
            )
    else:
        port = find_available_port(host)

    # Setup app
    cellxgene_url = f"http://{host}:{port}"

    # Import Flask app
    server = Server()

    server.create_app()
    server.app.config.update(SCRIPTS=scripts)

    if not verbose:
        log = logging.getLogger("werkzeug")
        log.setLevel(logging.ERROR)

    file_size = getsize(data)

    # if a big file, let the user know it may take a while to load.
    if file_size > BIG_FILE_SIZE_THRESHOLD:
        click.echo(f"[cellxgene] Loading data from {basename(data)}, this may take awhile...")
    else:
        click.echo(f"[cellxgene] Loading data from {basename(data)}.")

    # Fix for anaconda python. matplotlib typically expects python to be installed as a framework TKAgg is usually
    # available and fixes this issue. See https://matplotlib.org/faq/virtualenv_faq.html
    import matplotlib as mpl

    mpl.use("TkAgg")
    from server.app.scanpy_engine.scanpy_engine import ScanpyEngine

    try:
        server.attach_data(ScanpyEngine(data, e_args), title=title)
    except ScanpyFileError as e:
        raise click.ClickException(f"{e}")

    if open_browser:
        click.echo(f"[cellxgene] Launching! Opening your browser to {cellxgene_url} now.")
        webbrowser.open(cellxgene_url)
    else:
        click.echo(f"[cellxgene] Launching! Please go to {cellxgene_url} in your browser.")

    click.echo("[cellxgene] Type CTRL-C at any time to exit.")

    if not verbose:
        f = open(devnull, "w")
        sys.stdout = f

    try:
        server.app.run(host=host, debug=debug, port=port, threaded=True, use_debugger=False)
    except OSError as e:
        if e.errno == errno.EADDRINUSE:
            raise click.ClickException("Port is in use, please specify an open port using the --port flag.") from e
        raise
