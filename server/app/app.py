import os

from flask import Flask
from flask_caching import Cache
from flask_compress import Compress
from flask_cors import CORS

from server.app.rest_api.rest import get_api_resources
from server.app.util.utils import Float32JSONEncoder
from server.app.web import webapp

# Application Data
data = None
cache = Cache(config={"CACHE_TYPE": "simple", "CACHE_DEFAULT_TIMEOUT": 860_000})


def create_app():
    app = Flask(__name__, static_folder="web/static")
    app.json_encoder = Float32JSONEncoder
    cache.init_app(app)
    Compress(app)
    CORS(app)

    # Config
    SECRET_KEY = os.environ.get("CXG_SECRET_KEY", default="SparkleAndShine")
    app.config.update(SECRET_KEY=SECRET_KEY)


    resources = get_api_resources()
    app.register_blueprint(webapp.bp)
    app.register_blueprint(resources.blueprint)
    app.add_url_rule("/", endpoint="index")
    return app
