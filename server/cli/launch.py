import logging
from os import devnull
from os.path import splitext, basename
import sys
import warnings
import webbrowser

import click

from server.app.app import Server
from server.app.util.errors import ScanpyFileError
from server.app.util.utils import custom_format_warning
from server.utils.constants import MODES


@click.command()
@click.argument("data", metavar="<data file>", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option(
    "--layout",
    "-l",
    type=click.Choice(MODES),
    default="umap",
    show_default=True,
    help="Method for layout."
)
@click.option(
    "--diffexp",
    "-d",
    type=click.Choice(["ttest"]),
    default="ttest",
    show_default=True,
    help="Method for differential expression.",
)
@click.option("--title", "-t", help="Title to display (if omitted will use file name).", metavar="")
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
@click.option("--port", "-p", help="Port to run server on.", metavar="", default=5005, show_default=True)
@click.option("--obs-names", default=None, metavar="", help="Name of annotation field to use for observations.")
@click.option("--var-names", default=None, metavar="", help="Name of annotation to use for variables.")
@click.option("--host", default="127.0.0.1", help="Host IP address")
@click.option("--fixed-port", is_flag=True, default=False, show_default=True,
              help="Error if request port is not available instead of searching for another available port.")
@click.option(
    "--max-category-items",
    default=100,
    metavar="",
    show_default=True,
    help="Limits the number of categorical annotation items displayed.",
)
@click.option(
    "--diffexp-lfc-cutoff",
    default=0.01,
    show_default=True,
    help="Relative expression cutoff used when selecting top N differentially expressed genes",
)
@click.option(
    "--scripts",
    default=[],
    multiple=True,
    help="Additional script files to include in html page",
    show_default=True,
)
def launch(
        data,
        layout,
        diffexp,
        title,
        verbose,
        debug,
        obs_names,
        var_names,
        open_browser,
        port,
        host,
        fixed_port,
        max_category_items,
        diffexp_lfc_cutoff,
        scripts,
):
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

    if not verbose:
        sys.tracebacklimit = 0

    if not title:
        file_parts = splitext(basename(data))
        title = file_parts[0]

    if not fixed_port:
        new_port = find_available_port(host, port)
        if new_port != port:
            click.echo(
                f"[cellxgene] Warning: the port you specified was in use, using port {new_port} instead. "
                f"If you want to require cellxgene to only use the port specified and exit if that port is not "
                f"available run launch with the --fixed-port flag")
            port = new_port

    # Setup app
    cellxgene_url = f"http://{host}:{port}"

    # Import Flask app
    server = Server()

    server.create_app()
    server.app.config.update(SCRIPTS=scripts)

    if not verbose:
        log = logging.getLogger("werkzeug")
        log.setLevel(logging.ERROR)

    click.echo(f"[cellxgene] Loading data from {basename(data)}, this may take awhile...")

    # Fix for anaconda python. matplotlib typically expects python to be installed as a framework TKAgg is usually
    # available and fixes this issue. See https://matplotlib.org/faq/virtualenv_faq.html
    import matplotlib as mpl

    mpl.use("TkAgg")
    from server.app.scanpy_engine.scanpy_engine import ScanpyEngine

    args = {
        "layout": layout,
        "diffexp": diffexp,
        "max_category_items": max_category_items,
        "diffexp_lfc_cutoff": diffexp_lfc_cutoff,
        "obs_names": obs_names,
        "var_names": var_names,
    }

    try:
        server.attach_data(ScanpyEngine(data, args), title=title)
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
        server.app.run(host=host, debug=debug, port=port, threaded=True)
    except OSError as e:
        raise click.ClickException(f"{e}. Port is in use, please specify another port using the --port flag.")


def find_available_port(host, port):
    """
    Helper method to find open port on host. Tries 50 ports incremented from specified port
    """
    import errno
    import socket
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    new_port = port
    for i in range(3):
        new_port = port + i
        print(new_port)
        while True:
            try:
                s.bind((host, new_port))
                print("worked")
                break
            except socket.error as e:
                if e.errno == errno.EADDRINUSE:
                    continue
                else:
                    raise
            finally:
                s.close()
                break
    # else:
    #     click.echo(f"[cellxgene] No port in range {new_port} - {new_port + 50} was
    #     available to run cellxgene server on.")
    #     click.echo(f"[cellxgene] Exiting")
    #     raise OSError("No available ports")
    s.close()
    return new_port
