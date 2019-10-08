import errno
import functools
import logging
from os import devnull
from os.path import splitext, basename
import sys
import warnings
import webbrowser
from urllib.parse import urlparse

import click

from server.app.app import Server
from server.app.util.errors import ScanpyFileError
from server.app.util.utils import custom_format_warning
from server.utils.utils import find_available_port, is_port_available
from server.app.util.data_locator import DataLocator

# anything bigger than this will generate a special message
BIG_FILE_SIZE_THRESHOLD = 100 * 2 ** 20  # 100MB


def common_args(func):
    """
    Decorator to contain CLI args that will be common to both CLI and GUI: title and engine args.
    """

    @click.option("--title", "-t", help="Title to display (if omitted will use file name).")
    @click.option("--about",
                  help="A URL to more information about the dataset."
                       "(This must be an absolute URL including HTTP(S) protocol)")
    @click.option(
        "--embedding",
        "-e",
        default=[],
        multiple=True,
        show_default=False,
        help="Embedding name, eg, 'umap'. Repeat option for multiple embeddings. Defaults to all."
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
    @click.option(
        "--experimental-label-file",
        default=None,
        show_default=True,
        multiple=False,
        metavar="<user labels CSV file>",
        help="CSV file containing user annotations; will be overwritten.  Created if does not exist.",
    )
    @click.option(
        "--backed",
        is_flag=True,
        default=False,
        show_default=False,
        help="Load data in file-backed mode, which may save memory, but result in slower overall performance."
    )
    @click.option(
        "--disable-diffexp",
        is_flag=True,
        default=False,
        show_default=False,
        help="Disable on-demand differential expression."
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def parse_engine_args(embedding, obs_names, var_names, max_category_items,
                      diffexp_lfc_cutoff, experimental_label_file, backed,
                      disable_diffexp):
    return {
        "layout": embedding,
        "max_category_items": max_category_items,
        "diffexp_lfc_cutoff": diffexp_lfc_cutoff,
        "obs_names": obs_names,
        "var_names": var_names,
        "label_file": experimental_label_file,
        "backed": backed,
        "disable_diffexp": disable_diffexp
    }


@click.command()
@click.argument("data", nargs=1, metavar="<data file>", required=True)
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
        embedding,
        obs_names,
        var_names,
        max_category_items,
        diffexp_lfc_cutoff,
        title,
        scripts,
        about,
        experimental_label_file,
        backed,
        disable_diffexp
):
    """Launch the cellxgene data viewer.
    This web app lets you explore single-cell expression data.
    Data must be in a format that cellxgene expects, read the
    "getting started" guide.

    Examples:

    > cellxgene launch example_dataset/pbmc3k.h5ad --title pbmc3k

    > cellxgene launch <your data file> --title <your title>

    > cellxgene launch <url>"""

    e_args = parse_engine_args(embedding, obs_names, var_names, max_category_items,
                               diffexp_lfc_cutoff, experimental_label_file, backed,
                               disable_diffexp)
    try:
        data_locator = DataLocator(data)
    except RuntimeError as re:
        raise click.ClickException(f"Unable to access data at {data}.  {str(re)}")

    # Startup message
    click.echo("[cellxgene] Starting the CLI...")

    # Argument checking
    if data_locator.islocal():
        # if data locator is local, apply file system conventions and other "cheap"
        # validation checks.  If a URI, defer until we actually fetch the data and
        # try to read it.  Many of these tests don't make sense for URIs (eg, extension-
        # based typing).
        if not data_locator.exists():
            raise click.FileError(data, hint="file does not exist")
        if not data_locator.isfile():
            raise click.FileError(data, hint="data is not a file")
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

    if experimental_label_file:
        lf_name, lf_ext = splitext(experimental_label_file)
        if lf_ext and lf_ext != ".csv":
            raise click.FileError(basename(experimental_label_file), hint="label file type must be .csv")

    if about:
        def url_check(url):
            try:
                result = urlparse(url)
                if all([result.scheme, result.netloc]):
                    return True
                else:
                    return False
            except ValueError:
                return False

        if not url_check(about):
            raise click.ClickException("Must provide an absolute URL for --about. (Example format: http://example.com)")

    # Setup app
    cellxgene_url = f"http://{host}:{port}"

    # Import Flask app
    server = Server()

    server.create_app()
    server.app.config.update(SCRIPTS=scripts)

    if not verbose:
        log = logging.getLogger("werkzeug")
        log.setLevel(logging.ERROR)

    file_size = data_locator.size() if data_locator.islocal() else 0

    # if a big file, let the user know it may take a while to load.
    if file_size > BIG_FILE_SIZE_THRESHOLD:
        click.echo(f"[cellxgene] Loading data from {basename(data)}, this may take a while...")
    else:
        click.echo(f"[cellxgene] Loading data from {basename(data)}.")

    from server.app.scanpy_engine.scanpy_engine import ScanpyEngine

    try:
        server.attach_data(ScanpyEngine(data_locator, e_args), title=title, about=about)
    except ScanpyFileError as e:
        raise click.ClickException(f"{e}")

    if not disable_diffexp and server.app.data.config['diffexp_may_be_slow']:
        click.echo(f"[cellxgene] CAUTION: due to the size of your dataset, "
                   f"running differential expression may take longer or fail.")

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
        server.app.run(host=host, debug=debug, port=port, threaded=False if debug else True, use_debugger=False)
    except OSError as e:
        if e.errno == errno.EADDRINUSE:
            raise click.ClickException("Port is in use, please specify an open port using the --port flag.") from e
        raise
