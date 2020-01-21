import os
import datetime

from flask import Flask
from flask_caching import Cache
from flask_compress import Compress
from flask_cors import CORS

from server.app_single.rest import get_api_resources
from server.data_common.utils import Float32JSONEncoder
from server.common.web import webapp


class Server:
    def __init__(self, data, annotations, title="Demo", about=""):
        self.app = Flask(__name__, static_folder="../common/web/static")
        self.app.json_encoder = Float32JSONEncoder

        self.cache = Cache(config={"CACHE_TYPE": "simple", "CACHE_DEFAULT_TIMEOUT": 860_000})
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
        self.app.config.update(DATASET_TITLE=title, ABOUT_DATASET=about)

        self.app.data = data
        self.app.annotations = annotations
