import os

from flask import Flask
from flask_compress import Compress
from flask_cors import CORS
from flask_restful_swagger_2 import get_swagger_blueprint

from .web import webapp
from .rest_api.rest import get_api_resources

app = Flask(__name__)
Compress(app)
CORS(app)

# Config
CONFIG_FILE = os.environ.get("CXG_CONFIG_FILE", default="scanpy-test.cfg")
CXG_DIR = os.environ.get("CXG_DIRECTORY", default="/Users/charlotteweaver/Documents/Git/cxg-v2/data/")
SECRET_KEY = os.environ.get("CXG_SECRET_KEY", default="SparkleAndShine")
ENGINE = os.environ.get("CXG_ENGINE", default="scanpy")
TITLE = os.environ.get("DATASET_TITLE", default="PBMC 3K")
# TODO remove the 2 when this is prod
CXG_API_BASE = os.environ.get("CXG_API_BASE2", default="http://0.0.0.0:5005/api/")

if not CONFIG_FILE:
    raise ValueError("No config file set for Flask application")

# TODO check what is actually being configured here
app.config.from_pyfile(os.path.join(CXG_DIR, "config", CONFIG_FILE), silent=True)
app.config.update(
    SECRET_KEY=SECRET_KEY,
    CXG_API_BASE=CXG_API_BASE,
    ENGINE=ENGINE,
    DATA=CXG_DIR,
    DATASET_TITLE=TITLE
)

app.config['PROFILE'] = True
# app.wsgi_app = ProfilerMiddleware(app.wsgi_app, restrictions=[15])

# Application Data
data = None
if app.config["ENGINE"] == "scanpy":
    from .scanpy_engine.scanpy_engine import ScanpyEngine
    data = ScanpyEngine(app.config["DATA"], schema="data_schema.json")

REACTIVE_LIMIT = 1_000_000

# A list of swagger document objects
docs = []
resources = get_api_resources()
docs.append(resources.get_swagger_doc())


app.register_blueprint(webapp.bp)
app.register_blueprint(resources.blueprint)
app.register_blueprint(
    get_swagger_blueprint(docs, '/api/swagger', produces=["application/json"], title="cellxgene rest api",
                          description='An API connecting ExpressionMatrix2 clustering algorithm to cellxgene'))


app.add_url_rule('/', endpoint='index')
