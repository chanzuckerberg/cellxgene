import argparse
import logging
import os
import sys
import webbrowser

from flask import Flask
from flask_caching import Cache
from flask_compress import Compress
from flask_cors import CORS
from flask_restful_swagger_2 import get_swagger_blueprint

from .rest_api.rest import get_api_resources
from .util.utils import Float32JSONEncoder, whole_number
from .web import webapp

REACTIVE_LIMIT = 1_000_000

app = Flask(__name__, static_folder="web/static")
app.json_encoder = Float32JSONEncoder
cache = Cache(app, config={"CACHE_TYPE": "simple", "CACHE_DEFAULT_TIMEOUT": 860000})
Compress(app)
CORS(app)

# Config
SECRET_KEY = os.environ.get("CXG_SECRET_KEY", default="SparkleAndShine")

app.config.update(
    SECRET_KEY=SECRET_KEY,
)

# Application Data
data = None


# A list of swagger document objects
docs = []
resources = get_api_resources()
docs.append(resources.get_swagger_doc())

app.register_blueprint(webapp.bp)
app.register_blueprint(resources.blueprint)
app.register_blueprint(
    get_swagger_blueprint(docs, "/api/swagger", produces=["application/json"], title="cellxgene rest api",
                          description="An API connecting ExpressionMatrix2 clustering algorithm to cellxgene"))

app.add_url_rule("/", endpoint="index")


def create_cli():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.description = """
synopsis:
cellxgene <command> <data> [options]

description:

cellxgene is a local web application for exploring single cell expression.
    """
    parser.add_argument("-V", "--version", help="show version and exit")
    subparsers = parser.add_subparsers(dest="command")
    subparsers.required = True
    launch_group = subparsers.add_parser("launch", help="launch web application",
                                         formatter_class=argparse.RawTextHelpFormatter)
    launch_group.description = """
    cellxgene launches a local web application for exploring single cell expression data.

    Data must be in a format that cellxgene expects [[ how to format ]]
        """
    launch_group.epilog = """
    annotation names:
    The data viewer requires a unique, human readable name for each observation and variable.  These are used for
    various application features, such as the ability to view expression by gene. When launching cellxgene, appropriate
    observation and variable annotations must be identified.

    If --obs-name or --var-name parameters are specified, values in the named annotations will be used. If not
    specified, the observation and variable index values will name each respectively.  An error will generated if the
    values for each are not unique.

    examples:
    To run with the example dataset:

        cellxgene example_dataset/pbmc3k.h5ad --title PBMC3K

    To run with your own data with tsne layout:

        cellxgene <your data file> --title <your title> -l tsne

    To indicate that the human-readable variable annotation is named 'gene_names', and the human-readable observation
    is 'cell_names':

        cellxgene mydata.h5ad -var-name gene_names -obs-name cell_names

        """
    launch_group.add_argument("data", metavar="data", help="file containing the data to display")
    launch_group.add_argument("--title", "-t", help="title to display -- if this is omitted the title will be the name "
                                                    "of the data file.")
    launch_group.add_argument(
        "--listen-all",
        help="bind to all interfaces (this makes the server accessible beyond this computer)",
        action="store_true")
    launch_group.add_argument("--port", help="port to run server on", type=int, default=5005)
    launch_group.add_argument("-v", "--verbose", action="store_true",
                              help="more verbose output, including outputting warnings and every REST request")
    launch_group.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)
    launch_group.add_argument("--no-open", help="do not launch the webbrowser", action="store_false",
                              dest="open_browser")
    launch_group.add_argument(
        "--max-category-items",
        type=whole_number,
        help="maximum number of categories to display on the front-end. "
             "Annotations with more categories than this number will not be not displayed",
        default=100)
    try:
        from .scanpy_engine.scanpy_engine import ScanpyEngine
    except ImportError as e:
        # We will handle more engines when they come
        raise ImportError('Scanpy is required for cellxgene, please install scanpy and try again', e) from e
    else:
        ScanpyEngine.add_to_parser(launch_group)
    return parser


def run_scanpy(args):
    title = args.title
    if not title:
        file_parts = os.path.splitext(os.path.basename(args.data))
        title = file_parts[0]
    if args.listen_all:
        host = "0.0.0.0"
    else:
        host = "127.0.0.1"
    cellxgene_url = f"http://{host}:{args.port}"
    api_base = f"{cellxgene_url}/api/"
    app.config.update(
        DATASET_TITLE=title,
        CXG_API_BASE=api_base
    )

    if not args.verbose:
        log = logging.getLogger('werkzeug')
        log.setLevel(logging.ERROR)
    from .scanpy_engine.scanpy_engine import ScanpyEngine
    print(f"Loading data from {args.data} (this may take a while for large datasets)")
    app.data = ScanpyEngine(args.data, layout_method=args.layout, diffexp_method=args.diffexp,
                            max_category_items=args.max_category_items)
    print(f"Launching cellxgene")
    if args.open_browser:
        webbrowser.open(cellxgene_url)
    print(f"Please go to {cellxgene_url}")
    app.run(host=host, debug=args.debug, port=args.port)


def main():
    parser = create_cli()
    args = parser.parse_args()
    # Debug sets up developer mode
    if args.debug:
        args.verbose = True
        args.open_browser = False
    if not args.verbose:
        sys.tracebacklimit = 0
    # TODO pick engine based on input file
    print("cellxgene starting...\n")
    run_scanpy(args)
