import os
import datetime

from flask import Flask
from flask_caching import Cache
from flask_compress import Compress
from flask_cors import CORS

from server.app.rest_api.rest import get_api_resources
from server.app.util.utils import Float32JSONEncoder
from server.app.web import webapp


class Server:
    def __init__(self):
        self.data = None
        self.cache = Cache(config={"CACHE_TYPE": "simple", "CACHE_DEFAULT_TIMEOUT": 860_000})
        self.app = None

    def create_app(self):
        self.app = Flask(__name__, static_folder="web/static")
        self.app.json_encoder = Float32JSONEncoder
        self.cache.init_app(self.app)
        Compress(self.app)
        CORS(self.app, supports_credentials=True)

        # enable session data
        self.app.permanent_session_lifetime = datetime.timedelta(days=50 * 365)

        # Config
        SECRET_KEY = os.environ.get("CXG_SECRET_KEY", default="SparkleAndShine")
        self.app.config.update(SECRET_KEY=SECRET_KEY)
        self.app.config.update(SCRIPTS=[])

        resources = get_api_resources()
        self.app.register_blueprint(webapp.bp)
        self.app.register_blueprint(resources.blueprint)
        self.app.add_url_rule("/", endpoint="index")

    def attach_data(self, data, title="Demo", about=""):
        self.app.config.update(DATASET_TITLE=title, ABOUT_DATASET=about)
        self.app.data = data
