import argparse
import os

from flask import Flask
from flask_caching import Cache
from flask_compress import Compress
from flask_cors import CORS
from flask_restful_swagger_2 import get_swagger_blueprint

from .rest_api.rest import get_api_resources
from .web import webapp

REACTIVE_LIMIT = 1_000_000

app = Flask(__name__)
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
    global data
    title = args.title
    if not title:
        title = os.path.basename(os.path.normpath(args.data_directory))
    api_base = f"http://0.0.0.0:{args.port}/api/"
    app.config.update(
        DATASET_TITLE=title,
        CXG_API_BASE=api_base
    )

    from .scanpy_engine.scanpy_engine import ScanpyEngine
    data = ScanpyEngine(args.data_directory, schema="data_schema.json")
    app.run(host="0.0.0.0", debug=True, port=args.port)


def main():
    parser = argparse.ArgumentParser(description="Cellxgene is a tool for exploring single cell expression.")
    parser.add_argument("--title", "-t", help="Title to display -- if this is omitted the title will be the name "
                                              "of the directory from the data_directory arg")
    parser.add_argument("--port", help="Port to run server on.", type=int, default=5005)
    subparsers = parser.add_subparsers(dest="cellxgene_command")
    scanpy_subparser = subparsers.add_parser("scanpy", help="run cellxgene using the scanpy engine")
    scanpy_subparser.add_argument("data_directory", metavar="dir", help="Directory containing data and schema file")
    scanpy_subparser.set_defaults(func=run_scanpy)
    args = parser.parse_args()
    args.func(args)
