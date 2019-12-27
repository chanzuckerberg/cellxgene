import errno
import functools
import logging
from os import devnull, mkdir
from os.path import splitext, basename, isdir
import sys
import warnings
import webbrowser
from urllib.parse import urlparse

import click

from server.app.app import Server
from server.app.util.errors import ScanpyFileError
from server.app.util.utils import custom_format_warning
from server.utils.utils import find_available_port, is_port_available, sort_options
from server.app.util.data_locator import DataLocator

# anything bigger than this will generate a special message
BIG_FILE_SIZE_THRESHOLD = 100 * 2 ** 20  # 100MB


def common_args(func):
    """
    Decorator to contain CLI args that will be common to both CLI and GUI: title and engine args.
    """

    @click.option(
        "--title",
        "-t",
        metavar="<text>",
        help="Title to display. If omitted will use file name.")
    @click.option(
        "--about",
        metavar="<URL>",
        help="URL providing more information about the dataset "
             "(hint: must be a fully specified absolute URL).")
    @click.option(
        "--embedding",
        "-e",
        default=[],
        multiple=True,
        show_default=False,
        metavar="<text>",
        help="Embedding name, eg, 'umap'. Repeat option for multiple embeddings. Defaults to all."
    )
    @click.option(
        "--obs-names",
        "-obs",
        default=None,
        metavar="<text>",
        help="Name of annotation field to use for observations. If not specified cellxgene will use the the obs index.")
    @click.option(
        "--var-names",
        "-var",
        default=None,
        metavar="<text>",
        help="Name of annotation to use for variables. If not specified cellxgene will use the the var index.")
    @click.option(
        "--max-category-items",
        default=1000,
        metavar="<integer>",
        show_default=True,
        help="Will not display categories with more distinct values than specified.",)
    @click.option(
        "--diffexp-lfc-cutoff",
        "-de",
        default=0.01,
        show_default=True,
        metavar="<float>",
        help="Minimum log fold change threshold for differential expression.",)
    @click.option(
        "--experimental-annotations",
        is_flag=True,
        default=False,
        show_default=True,
        help="Enable user annotation of data."
    )
    @click.option(
        "--experimental-annotations-file",
        default=None,
        show_default=True,
        multiple=False,
        metavar="<path>",
        help="CSV file to initialize editing of existing annotations; will be altered in-place. "
             "Incompatible with --annotations-output-dir.",)
    @click.option(
        "--experimental-annotations-output-dir",
        default=None,
        show_default=False,
        multiple=False,
        metavar="<directory path>",
        help="Directory of where to save output annotations; filename will be specified in the application. "
             "Incompatible with --annotations-input-file.",)
    @click.option(
        "--backed",
        "-b",
        is_flag=True,
        default=False,
        show_default=False,
        help="Load data in file-backed mode. This may save memory, but may result in slower overall performance.")
    @click.option(
        "--disable-diffexp",
        is_flag=True,
        default=False,
        show_default=False,
        help="Disable on-demand differential expression.")
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def parse_engine_args(embedding, obs_names, var_names, max_category_items, diffexp_lfc_cutoff,
                      experimental_annotations, experimental_annotations_file,
                      experimental_annotations_output_dir, backed, disable_diffexp):
    annotations_file = experimental_annotations_file if experimental_annotations else None
    annotations_output_dir = experimental_annotations_output_dir if experimental_annotations else None
    return {
        "layout": embedding,
        "max_category_items": max_category_items,
        "diffexp_lfc_cutoff": diffexp_lfc_cutoff,
        "obs_names": obs_names,
        "var_names": var_names,
        "annotations": experimental_annotations,
        "annotations_file": annotations_file,
        "annotations_output_dir": annotations_output_dir,
        "backed": backed,
        "disable_diffexp": disable_diffexp
    }


@sort_options
@click.command(short_help="Launch the cellxgene data viewer. "
                          "Run `cellxgene launch --help` for more information.",
               options_metavar="<options>",)
@click.argument("data", nargs=1, metavar="<path to data file>", required=True)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    show_default=True,
    help="Provide verbose output, including warnings and all server requests.",)
@click.option(
    "--debug",
    "-d",
    is_flag=True,
    default=False,
    show_default=True,
    help="Run in debug mode. This is helpful for cellxgene developers, "
         "or when you want more information about an error condition.",)
@click.option(
    "--open",
    "-o",
    "open_browser",
    is_flag=True,
    default=False,
    show_default=True,
    help="Open web browser after launch.",)
@click.option(
    "--port",
    "-p",
    metavar="<port>",
    show_default=True,
    help="Port to run server on. If not specified cellxgene will find an available port.",)
@click.option(
    "--host",
    metavar="<IP address>",
    default="127.0.0.1",
    show_default=False,
    help="Host IP address. By default cellxgene will use localhost (e.g. 127.0.0.1).")
@click.option(
    "--scripts",
    "-s",
    default=[],
    multiple=True,
    metavar="<text>",
    help="Additional script files to include in HTML page. If not specified, "
         "no additional script files will be included.",
    show_default=False,)
@click.help_option("--help", "-h", help="Show this message and exit.")
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
        experimental_annotations,
        experimental_annotations_file,
        experimental_annotations_output_dir,
        backed,
        disable_diffexp
):
    """Launch the cellxgene data viewer.
    This web app lets you explore single-cell expression data.
    Data must be in a format that cellxgene expects.
    Read the "getting started" guide to learn more:
    https://chanzuckerberg.github.io/cellxgene/getting-started.html

    Examples:

    > cellxgene launch example_dataset/pbmc3k.h5ad --title pbmc3k

    > cellxgene launch <your data file> --title <your title>

    > cellxgene launch <url>"""

    e_args = parse_engine_args(embedding, obs_names, var_names, max_category_items,
                               diffexp_lfc_cutoff,
                               experimental_annotations,
                               experimental_annotations_file,
                               experimental_annotations_output_dir,
                               backed,
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

    if not experimental_annotations:
        if experimental_annotations_file is not None:
            click.echo("Warning: --experimental-annotations-file ignored as --annotations not enabled.")
        if experimental_annotations_output_dir is not None:
            click.echo("Warning: --experimental-annotations-output-dir ignored as --annotations not enabled.")
    else:
        if experimental_annotations_file is not None and experimental_annotations_output_dir is not None:
            raise click.ClickException("--experimental-annotations-file and --experimental-annotations-output-dir "
                                       "may not be used together.")

        if experimental_annotations_file is not None:
            lf_name, lf_ext = splitext(experimental_annotations_file)
            if lf_ext and lf_ext != ".csv":
                raise click.FileError(basename(experimental_annotations_file), hint="annotation file type must be .csv")

        if experimental_annotations_output_dir is not None and not isdir(experimental_annotations_output_dir):
            try:
                mkdir(experimental_annotations_output_dir)
            except OSError:
                raise click.ClickException("Unable to create directory specified by "
                                           "--experimental-annotations-output-dir")

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
