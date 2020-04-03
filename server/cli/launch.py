import errno
import functools
import logging
from os import devnull
import sys
import webbrowser

import click
from flask_compress import Compress
from flask_cors import CORS

from server.common.utils import sort_options
from server.common.errors import DatasetAccessError, ConfigurationError
from server.common.app_config import AppConfig
from server.common.default_config import default_config
from server.app.app import Server

DEFAULT_CONFIG = AppConfig()


def annotation_args(func):
    @click.option(
        "--disable-annotations",
        is_flag=True,
        default=not DEFAULT_CONFIG.user_annotations__enable,
        show_default=True,
        help="Disable user annotation of data.",
    )
    @click.option(
        "--annotations-file",
        default=DEFAULT_CONFIG.user_annotations__local_file_csv__file,
        show_default=True,
        multiple=False,
        metavar="<path>",
        help="CSV file to initialize editing of existing annotations; will be altered in-place. "
        "Incompatible with --annotations-dir.",
    )
    @click.option(
        "--annotations-dir",
        default=DEFAULT_CONFIG.user_annotations__local_file_csv__directory,
        show_default=False,
        multiple=False,
        metavar="<directory path>",
        help="Directory of where to save output annotations; filename will be specified in the application. "
        "Incompatible with --annotations-file.",
    )
    @click.option(
        "--experimental-annotations-ontology",
        is_flag=True,
        default=DEFAULT_CONFIG.user_annotations__ontology__enable,
        show_default=True,
        help="When creating annotations, optionally autocomplete names from ontology terms.",
    )
    @click.option(
        "--experimental-annotations-ontology-obo",
        default=DEFAULT_CONFIG.user_annotations__ontology__obo_location,
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
        default=DEFAULT_CONFIG.presentation__max_categories,
        metavar="<integer>",
        show_default=True,
        help="Will not display categories with more distinct values than specified.",
    )
    @click.option(
        "--diffexp-lfc-cutoff",
        "-de",
        default=DEFAULT_CONFIG.diffexp__lfc_cutoff,
        show_default=True,
        metavar="<float>",
        help="Minimum log fold change threshold for differential expression.",
    )
    @click.option(
        "--disable-diffexp",
        is_flag=True,
        default=not DEFAULT_CONFIG.diffexp__enable,
        show_default=False,
        help="Disable on-demand differential expression.",
    )
    @click.option(
        "--embedding",
        "-e",
        default=DEFAULT_CONFIG.embeddings__names,
        multiple=True,
        show_default=False,
        metavar="<text>",
        help="Embedding name, eg, 'umap'. Repeat option for multiple embeddings. Defaults to all.",
    )
    @click.option(
        "--experimental-enable-reembedding",
        is_flag=True,
        default=DEFAULT_CONFIG.embeddings__enable_reembedding,
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
        default=DEFAULT_CONFIG.single_dataset__obs_names,
        metavar="<text>",
        help="Name of annotation field to use for observations. If not specified cellxgene will use the the obs index.",
    )
    @click.option(
        "--var-names",
        "-var",
        default=DEFAULT_CONFIG.single_dataset__var_names,
        metavar="<text>",
        help="Name of annotation to use for variables. If not specified cellxgene will use the the var index.",
    )
    @click.option(
        "--backed",
        "-b",
        is_flag=True,
        default=DEFAULT_CONFIG.adaptor__anndata_adaptor__backed,
        show_default=False,
        help="Load anndata in file-backed mode. " "This may save memory, but may result in slower overall performance.",
    )
    @click.option(
        "--title",
        "-t",
        default=DEFAULT_CONFIG.single_dataset__title,
        metavar="<text>",
        help="Title to display. If omitted will use file name.",
    )
    @click.option(
        "--about",
        default=DEFAULT_CONFIG.single_dataset__about,
        metavar="<URL>",
        help="URL providing more information about the dataset (hint: must be a fully specified absolute URL).",
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
        default=DEFAULT_CONFIG.server__debug,
        show_default=True,
        help="Run in debug mode. This is helpful for cellxgene developers, "
        "or when you want more information about an error condition.",
    )
    @click.option(
        "--verbose",
        "-v",
        is_flag=True,
        default=DEFAULT_CONFIG.server__verbose,
        show_default=True,
        help="Provide verbose output, including warnings and all server requests.",
    )
    @click.option(
        "--port",
        "-p",
        metavar="<port>",
        default=DEFAULT_CONFIG.server__port,
        type=int,
        show_default=True,
        help="Port to run server on. If not specified cellxgene will find an available port.",
    )
    @click.option(
        "--host",
        metavar="<IP address>",
        default=DEFAULT_CONFIG.server__host,
        show_default=False,
        help="Host IP address. By default cellxgene will use localhost (e.g. 127.0.0.1).",
    )
    @click.option(
        "--scripts",
        "-s",
        default=DEFAULT_CONFIG.server__scripts,
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
        default=DEFAULT_CONFIG.multi_dataset__dataroot,
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
        default=DEFAULT_CONFIG.server__open_browser,
        show_default=True,
        help="Open web browser after launch.",
    )
    @click.option(
        "--config-file",
        "-c",
        "config_file",
        default=None,
        show_default=True,
        help="Location to yaml file with configuration settings",
    )
    @click.option(
        "--dump-default-config",
        "dump_default_config",
        is_flag=True,
        default=False,
        show_default=True,
        help="Print default configuration settings and exit",
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


class CliLaunchServer(Server):
    """
    the CLI runs a local web server, and needs to enable a few more features.
    """

    def __init__(self, app_config):
        super().__init__(app_config)

    @staticmethod
    def _before_adding_routes(app, app_config):
        app.config["COMPRESS_MIMETYPES"] = [
            "text/html",
            "text/css",
            "text/xml",
            "application/json",
            "application/javascript",
            "application/octet-stream",
        ]
        Compress(app)
        if app_config.server__debug:
            CORS(app, supports_credentials=True)


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
    config_file,
    dump_default_config,
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

    if dump_default_config:
        print(default_config)
        sys.exit(0)

    # Startup message
    click.echo("[cellxgene] Starting the CLI...")

    # app config
    app_config = AppConfig()

    try:
        if config_file:
            app_config.update_from_config_file(config_file)

        # Determine which config options were give on the command line.
        # Those will override the ones provided in the config file (if provided).
        cli_config = AppConfig()
        cli_config.update(
            server__verbose=verbose,
            server__debug=debug,
            server__host=host,
            server__port=port,
            server__scripts=scripts,
            server__open_browser=open_browser,
            single_dataset__datapath=datapath,
            single_dataset__title=title,
            single_dataset__about=about,
            single_dataset__obs_names=obs_names,
            single_dataset__var_names=var_names,
            multi_dataset__dataroot=dataroot,
            user_annotations__enable=not disable_annotations,
            user_annotations__local_file_csv__file=annotations_file,
            user_annotations__local_file_csv__directory=annotations_dir,
            user_annotations__ontology__enable=experimental_annotations_ontology,
            user_annotations__ontology__obo_location=experimental_annotations_ontology_obo,
            presentation__max_categories=max_category_items,
            embeddings__names=embedding,
            embeddings__enable_reembedding=experimental_enable_reembedding,
            diffexp__enable=not disable_diffexp,
            diffexp__lfc_cutoff=diffexp_lfc_cutoff,
            adaptor__anndata_adaptor__backed=backed,
        )
        diff = cli_config.changes_from_default()
        changes = {}
        for key, val, defval in diff:
            changes[key] = val
        app_config.update(**changes)

        # process the configuration
        #  any errors will be thrown as an exception.
        #  any info messages will be passed to the messagefn function.

        def messagefn(message):
            click.echo("[cellxgene] " + message)

        # Use a default secret if one is not provided
        if not app_config.server__flask_secret_key:
            app_config.update(server__flask_secret_key="SparkleAndShine")

        app_config.complete_config(messagefn)

    except (ConfigurationError, DatasetAccessError) as e:
        raise click.ClickException(e)

    handle_scripts(scripts)

    # create the server
    server = CliLaunchServer(app_config)

    if not app_config.server__verbose:
        log = logging.getLogger("werkzeug")
        log.setLevel(logging.ERROR)

    cellxgene_url = f"http://{app_config.server__host}:{app_config.server__port}"
    if app_config.server__open_browser:
        click.echo(f"[cellxgene] Launching! Opening your browser to {cellxgene_url} now.")
        webbrowser.open(cellxgene_url)
    else:
        click.echo(f"[cellxgene] Launching! Please go to {cellxgene_url} in your browser.")

    click.echo("[cellxgene] Type CTRL-C at any time to exit.")

    if not app_config.server__verbose:
        f = open(devnull, "w")
        sys.stdout = f

    try:
        server.app.run(
            host=app_config.server__host,
            debug=app_config.server__debug,
            port=app_config.server__port,
            threaded=not app_config.server__debug,
            use_debugger=False,
            use_reloader=False,
        )
    except OSError as e:
        if e.errno == errno.EADDRINUSE:
            raise click.ClickException("Port is in use, please specify an open port using the --port flag.") from e
        raise
