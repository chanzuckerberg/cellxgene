import errno
import functools
import logging
from os import devnull, mkdir, environ
from os.path import splitext, basename, isdir
import sys
import warnings
import webbrowser
from urllib.parse import urlparse

import click

from server.common.utils import custom_format_warning
from server.common.utils import find_available_port, is_port_available, sort_options
from server.common.errors import DatasetAccessError
from server.data_common.matrix_loader import MatrixDataLoader, MatrixDataCacheManager, MatrixDataType
from server.common.annotations import AnnotationsLocalFile
from server.common.app_config import AppConfig

from server.common.errors import OntologyLoadFailure

# anything bigger than this will generate a special message
BIG_FILE_SIZE_THRESHOLD = 100 * 2 ** 20  # 100MB
DEFAULT_SERVER_PORT = int(environ.get("CXG_SERVER_PORT", "5005"))


def annotation_args(func):
    @click.option(
        "--disable-annotations",
        is_flag=True,
        default=False,
        show_default=True,
        help="Disable user annotation of data.",
    )
    @click.option(
        "--annotations-file",
        default=None,
        show_default=True,
        multiple=False,
        metavar="<path>",
        help="CSV file to initialize editing of existing annotations; will be altered in-place. "
        "Incompatible with --annotations-dir.",
    )
    @click.option(
        "--annotations-dir",
        default=None,
        show_default=False,
        multiple=False,
        metavar="<directory path>",
        help="Directory of where to save output annotations; filename will be specified in the application. "
        "Incompatible with --annotations-file.",
    )
    @click.option(
        "--experimental-annotations-ontology",
        is_flag=True,
        default=False,
        show_default=True,
        help="When creating annotations, optionally autocomplete names from ontology terms.",
    )
    @click.option(
        "--experimental-annotations-ontology-obo",
        default=None,
        show_default=True,
        metavar="<path or url>",
        help="Location of OBO file defining cell annotation autosuggest terms.",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def config_args(func):
    @click.option(
        "--max-category-items",
        default=1000,
        metavar="<integer>",
        show_default=True,
        help="Will not display categories with more distinct values than specified.",
    )
    @click.option(
        "--diffexp-lfc-cutoff",
        "-de",
        default=0.01,
        show_default=True,
        metavar="<float>",
        help="Minimum log fold change threshold for differential expression.",
    )
    @click.option(
        "--disable-diffexp",
        is_flag=True,
        default=False,
        show_default=False,
        help="Disable on-demand differential expression.",
    )
    @click.option(
        "--embedding",
        "-e",
        default=[],
        multiple=True,
        show_default=False,
        metavar="<text>",
        help="Embedding name, eg, 'umap'. Repeat option for multiple embeddings. Defaults to all.",
    )
    @click.option(
        "--experimental-enable-reembedding",
        is_flag=True,
        default=False,
        show_default=False,
        hidden=True,
        help="Enable experimental on-demand re-embedding using UMAP. WARNING: may be very slow.",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def dataset_args(func):
    @click.option(
        "--obs-names",
        "-obs",
        default=None,
        metavar="<text>",
        help="Name of annotation field to use for observations. If not specified cellxgene will use the the obs index.",
    )
    @click.option(
        "--var-names",
        "-var",
        default=None,
        metavar="<text>",
        help="Name of annotation to use for variables. If not specified cellxgene will use the the var index.",
    )
    @click.option(
        "--backed",
        "-b",
        is_flag=True,
        default=False,
        show_default=False,
        help="Load anndata in file-backed mode. " "This may save memory, but may result in slower overall performance.",
    )
    @click.option("--title", "-t", metavar="<text>", help="Title to display. If omitted will use file name.")
    @click.option(
        "--about",
        metavar="<URL>",
        help="URL providing more information about the dataset " "(hint: must be a fully specified absolute URL).",
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def server_args(func):
    @click.option(
        "--debug",
        "-d",
        is_flag=True,
        default=False,
        show_default=True,
        help="Run in debug mode. This is helpful for cellxgene developers, "
        "or when you want more information about an error condition.",
    )
    @click.option(
        "--verbose",
        "-v",
        is_flag=True,
        default=False,
        show_default=True,
        help="Provide verbose output, including warnings and all server requests.",
    )
    @click.option(
        "--port",
        "-p",
        metavar="<port>",
        default=None,
        show_default=True,
        help="Port to run server on. If not specified cellxgene will find an available port.",
    )
    @click.option(
        "--host",
        metavar="<IP address>",
        default="127.0.0.1",
        show_default=False,
        help="Host IP address. By default cellxgene will use localhost (e.g. 127.0.0.1).",
    )
    @click.option(
        "--scripts",
        "-s",
        default=[],
        multiple=True,
        metavar="<text>",
        help="Additional script files to include in HTML page. If not specified, "
        "no additional script files will be included.",
        show_default=False,
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def launch_args(func):
    @annotation_args
    @config_args
    @dataset_args
    @server_args
    @click.option(
        "--dataroot",
        default=None,
        metavar="<data directory>",
        help="Enable cellxgene to serve multiple files. Supply path (local directory or URL)"
        " to folder containing H5AD and/or CXG datasets.",
        hidden=True,
    )  # TODO, unhide when dataroot is supported)
    @click.argument("datapath", required=False, metavar="<path to data file>")
    @click.option(
        "--open",
        "-o",
        "open_browser",
        is_flag=True,
        default=False,
        show_default=True,
        help="Open web browser after launch.",
    )
    @click.help_option("--help", "-h", help="Show this message and exit.")
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def handle_scripts(scripts):
    if scripts:
        click.echo(
            r"""
    / / /\ \ \__ _ _ __ _ __ (_)_ __   __ _
    \ \/  \/ / _` | '__| '_ \| | '_ \ / _` |
     \  /\  / (_| | |  | | | | | | | | (_| |
      \/  \/ \__,_|_|  |_| |_|_|_| |_|\__, |
                                      |___/
    The --scripts flag is intended for developers to include google analytics etc. You could be opening yourself to a
    security risk by including the --scripts flag. Make sure you trust the scripts that you are including.
            """
        )
        scripts_pretty = ", ".join(scripts)
        click.confirm(f"Are you sure you want to inject these scripts: {scripts_pretty}?", abort=True)


def handle_verbose(verbose):
    if not verbose:
        sys.tracebacklimit = 0


@sort_options
@click.command(
    short_help="Launch the cellxgene data viewer. " "Run `cellxgene launch --help` for more information.",
    options_metavar="<options>",
)
@launch_args
def launch(
    datapath,
    dataroot,
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
    disable_annotations,
    annotations_file,
    annotations_dir,
    backed,
    disable_diffexp,
    experimental_annotations_ontology,
    experimental_annotations_ontology_obo,
    experimental_enable_reembedding,
):
    """Launch the cellxgene data viewer.
    This web app lets you explore single-cell expression data.
    Data must be in a format that cellxgene expects.
    Read the "getting started" guide to learn more:
    https://chanzuckerberg.github.io/cellxgene/getting-started.html

    Examples:

    > cellxgene launch example-dataset/pbmc3k.h5ad --title pbmc3k

    > cellxgene launch <your data file> --title <your title>

    > cellxgene launch <url>"""

    # TODO Examples to provide when "--dataroot" is unhidden
    # > cellxgene launch --dataroot example-dataset/
    #
    # > cellxgene launch --dataroot <url>

    # Startup message
    click.echo("[cellxgene] Starting the CLI...")

    if datapath is None and dataroot is None:
        # TODO:  change the error message once dataroot is fully supported
        raise click.ClickException('Missing argument "<path to data file>."')
        # raise click.ClickException("must supply either <path to data file> or --dataroot")
    if datapath is not None and dataroot is not None:
        raise click.ClickException("must supply only one of <path to data file> or --dataroot")

    if datapath:
        # preload this data set
        matrix_data_loader = MatrixDataLoader(datapath)

        try:
            matrix_data_loader.pre_load_validation()
        except DatasetAccessError as e:
            raise click.ClickException(str(e))

        if experimental_enable_reembedding:
            if matrix_data_loader.matrix_data_type() != MatrixDataType.H5AD:
                raise click.ClickException("--experimental-enable-reembedding is only supported with H5AD files.")
            if backed:
                raise click.ClickException(
                    "--experimental-enable-reembedding is not supported when run in --backed mode."
                )

        file_size = matrix_data_loader.file_size()
        if file_size > BIG_FILE_SIZE_THRESHOLD:
            click.echo(f"[cellxgene] Loading data from {basename(datapath)}, this may take a while...")
        else:
            click.echo(f"[cellxgene] Loading data from {basename(datapath)}.")

    if debug:
        verbose = True
        open_browser = False
    else:
        warnings.formatwarning = custom_format_warning

    handle_verbose(verbose)
    handle_scripts(scripts)

    if port:
        if debug:
            raise click.ClickException("--port and --debug may not be used together (try --verbose for error logging).")
        if not is_port_available(host, int(port)):
            raise click.ClickException(
                f"The port selected {port} is in use, please specify an open port using the --port flag."
            )
    else:
        port = find_available_port(host, DEFAULT_SERVER_PORT)

    if disable_annotations:
        if annotations_file is not None:
            click.echo("Warning: --annotations-file ignored as annotations are disabled.")
        if annotations_dir is not None:
            click.echo("Warning: --annotations-dir ignored as annotations are disabled.")
        if experimental_annotations_ontology:
            click.echo("Warning: --experimental-annotations-ontology ignored as annotations are disabled.")
        if experimental_annotations_ontology_obo is not None:
            click.echo("Warning: --experimental-annotations-ontology-obo ignored as annotations are disabled.")
    else:
        if annotations_file is not None and annotations_dir is not None:
            raise click.ClickException(
                "--annotations-file and --annotations-dir " "may not be used together."
            )

        if annotations_file is not None:
            lf_name, lf_ext = splitext(annotations_file)
            if lf_ext and lf_ext != ".csv":
                raise click.FileError(basename(annotations_file), hint="annotation file type must be .csv")

        if annotations_dir is not None and not isdir(annotations_dir):
            try:
                mkdir(annotations_dir)
            except OSError:
                raise click.ClickException(
                    "Unable to create directory specified by " "--annotations-dir"
                )

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

    # app config
    app_config = AppConfig(
        datapath=datapath,
        dataroot=dataroot,
        title=title,
        about=about,
        scripts=scripts,
        layout=embedding,
        max_category_items=max_category_items,
        diffexp_lfc_cutoff=diffexp_lfc_cutoff,
        obs_names=obs_names,
        var_names=var_names,
        anndata_backed=backed,
        disable_diffexp=disable_diffexp,
        enable_reembedding=experimental_enable_reembedding,
    )

    matrix_data_cache_manager = MatrixDataCacheManager()
    data_adaptor = None
    if datapath:
        try:
            with matrix_data_cache_manager.data_adaptor(datapath, app_config) as data_adaptor:
                if not disable_diffexp and data_adaptor.parameters.get("diffexp_may_be_slow", False):
                    click.echo(
                        f"[cellxgene] CAUTION: due to the size of your dataset, "
                        f"running differential expression may take longer or fail."
                    )
        except Exception as e:
            raise click.ClickException(str(e))

    # create an annotations object.  Only AnnotationsLocalFile is used (for now)
    annotations = None

    if not disable_annotations:
        annotations = AnnotationsLocalFile(annotations_dir, annotations_file)

        # if the user has specified a fixed label file, go ahead and validate it
        # so that we can remove errors early in the process.

        if annotations_file and data_adaptor:
            data_adaptor.check_new_labels(annotations.read_labels(data_adaptor))

        if experimental_annotations_ontology or bool(experimental_annotations_ontology_obo):
            try:
                annotations.load_ontology(experimental_annotations_ontology_obo)
            except OntologyLoadFailure as e:
                raise click.ClickException("Unable to load ontology terms\n" + str(e))

    # create the server
    from server.app.app import Server

    server = Server(matrix_data_cache_manager, annotations, app_config)

    if not verbose:
        log = logging.getLogger("werkzeug")
        log.setLevel(logging.ERROR)

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
