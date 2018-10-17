import argparse
import os
import warnings
import webbrowser

from flask import Flask
from flask_caching import Cache
from flask_compress import Compress
from flask_cors import CORS
from flask_restful_swagger_2 import get_swagger_blueprint

from .rest_api.rest import get_api_resources
from .util.utils import Float32JSONEncoder
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


def run_scanpy(args):
    title = args.title
    if not title:
        title = os.path.basename(os.path.normpath(args.data_directory))
    api_base = f"http://127.0.0.1:{args.port}/api/"
    app.config.update(
        DATASET_TITLE=title,
        CXG_API_BASE=api_base
    )

    from .scanpy_engine.scanpy_engine import ScanpyEngine
    app.data = ScanpyEngine(args.data_directory, layout_method=args.layout, diffexp_method=args.diffexp)
    if args.bind_all:
        host = "0.0.0.0"
    else:
        host = "127.0.0.1"
    if args.open_browser:
        webbrowser.open(f"http://{host}:{args.port}")
    app.run(host=host, debug=True, port=args.port)


def main():
    parser = argparse.ArgumentParser(description="Cellxgene is a tool for exploring single cell expression.")
    parser.add_argument("--title", "-t", help="Title to display -- if this is omitted the title will be the name "
                                              "of the directory from the data_directory arg")
    parser.add_argument("--port", help="Port to run server on.", type=int, default=5005)
    parser.add_argument("--no-open", help="Do not launch the webbrowser.", action="store_false", dest="open_browser")
    parser.add_argument(
        "--bind-all",
        help="Bind to all interfaces (this makes the server accessible beyond this computer)",
        action="store_true")
    subparsers = parser.add_subparsers(dest="engine")
    subparsers.required = True
    try:
        from .scanpy_engine.scanpy_engine import ScanpyEngine
    except ImportError:
        warnings.simplefilter('default', ImportWarning)  # Enable ImportWarning
        warnings.warn("Scanpy engine not available", ImportWarning)
    else:
        ScanpyEngine.add_to_parser(subparsers, run_scanpy)
    if len(subparsers.choices) == 0:
        raise ImportError('Could not import any engines, see warnings above')
    args = parser.parse_args()
    args.func(args)
