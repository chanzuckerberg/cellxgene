import os
import argparse

from flask import Flask
from flask_compress import Compress
from flask_cors import CORS
from flask_restful_swagger_2 import get_swagger_blueprint

from .web import webapp
from .rest_api.rest import get_api_resources

REACTIVE_LIMIT = 1_000_000

app = Flask(__name__)
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


def run(args):
    global data
    title = args.title
    if not title:
        title = os.path.basename(os.path.normpath(args.data_directory))
    api_base = f"http://0.0.0.0:{args.port}/api/"
    app.config.update(
        DATASET_TITLE=title,
        CXG_API_BASE=api_base
    )
    if args.engine == "scanpy":
        from .scanpy_engine.scanpy_engine import ScanpyEngine
        data = ScanpyEngine(args.data_directory, schema="data_schema.json")
    app.run(host="0.0.0.0", debug=True, port=args.port)


def main():
    parser = argparse.ArgumentParser(description="Cellxgene is a tool for exploring single cell expression.")
    subparsers = parser.add_subparsers(dest="cellxgene_command")
    run_subparser = subparsers.add_parser("run", help="run cellxgene")
    run_subparser.add_argument("--title", "-t", help="Title to display (default = data directory name")
    run_subparser.add_argument("--port", help="Port to run server on.", type=int, default=5005)
    run_subparser.add_argument("engine", help="The underlying structure of the data")
    run_subparser.add_argument("data_directory", metavar="dir", help="Directory containing data and schema")
    run_subparser.set_defaults(func=run)
    args = parser.parse_args()

    args.func(args)

    # app.run(host="0.0.0.0", debug=True, port=5005)
