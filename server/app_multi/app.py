import os
import datetime

from flask import Flask, redirect, current_app, make_response, render_template
from flask import Blueprint, request
from flask_caching import Cache
from flask_compress import Compress
from flask_cors import CORS
from flask_restful import Api, Resource

from server.data_common.utils import Float32JSONEncoder, MatrixDataLoader

from http import HTTPStatus
from functools import wraps
import server.common.rest as common_rest
from server.common.utils import path_join

web_bp = Blueprint("webapp", "server.common.web", template_folder="templates")


@web_bp.route("/")
def index():
    # FIXME with a splash screen that includes a listing of all the datasets.
    # or perhaps a login screen if this is a hosted environment
    return "<H1>Welcome to cellxgene<H1>"


def get_data_engine(dataset):
    config = current_app.app_config
    location = path_join(config.dataroot, dataset)
    matrix_data_loader = MatrixDataLoader(location)
    matrix_data_loader.pre_checks()
    data = matrix_data_loader.open(current_app.app_config)
    return data


def rest_get_data_engine(func):
    @wraps(func)
    def wrapped_function(self, dataset):
        try:
            data_engine = get_data_engine(dataset)
            return func(self, data_engine)
        except RuntimeError as e:
            return make_response(f"Invalid dataset {dataset}: {str(e)}", HTTPStatus.BAD_REQUEST)
    return wrapped_function


def dataset_index(dataset):
    config = current_app.app_config
    location = path_join(config.dataroot, dataset)
    matrix_data_loader = MatrixDataLoader(location)
    try:
        matrix_data_loader.pre_checks()
        return render_template("index.html")
    except RuntimeError as e:
        return make_response(f"CXG load error {dataset}: {str(e)}", HTTPStatus.BAD_REQUEST)
    except Exception as e:
        return make_response(f"CXG load error {dataset}: {str(e)}", HTTPStatus.BAD_REQUEST)


def static_redirect(dataset, therest):
    """ redirect all static requests to the standard location """
    return redirect(f'/static/{therest}', code=301)


def favicon_redirect(dataset):
    """ redirect favicon to static dir """
    return redirect('/static/favicon.png', code=301)


class SchemaAPI(Resource):
    @rest_get_data_engine
    def get(self, data_engine):
        return common_rest.schema_get(data_engine, current_app.annotations)


class ConfigAPI(Resource):
    @rest_get_data_engine
    def get(self, data_engine):
        return common_rest.config_get(
            current_app.app_config, data_engine, current_app.annotations)


class AnnotationsObsAPI(Resource):
    @rest_get_data_engine
    def get(self, data_engine):
        return common_rest.annotations_obs_get(
            request, data_engine, current_app.annotations)

    @rest_get_data_engine
    def put(self, data_engine):
        return common_rest.annotations_obs_put(
            request, data_engine, current_app.annotations)


class AnnotationsVarAPI(Resource):
    @rest_get_data_engine
    def get(self, data_engine):
        return common_rest.annotations_var_get(request, data_engine, current_app.annotations)


class DataVarAPI(Resource):
    @rest_get_data_engine
    def put(self, data_engine):
        return common_rest.data_var_put(request, data_engine)


class DiffExpObsAPI(Resource):
    @rest_get_data_engine
    def post(self, data_engine):
        return common_rest.diffexp_obs_post(request, data_engine)


class LayoutObsAPI(Resource):
    @rest_get_data_engine
    def get(self, data_engine):
        return common_rest.layout_obs_get(request, data_engine)


def get_api_resources():
    bp = Blueprint("api", __name__, url_prefix="/<dataset>/api/v0.2")
    api = Api(bp)
    # Initialization routes
    api.add_resource(SchemaAPI, "/schema")
    api.add_resource(ConfigAPI, "/config")
    # Data routes
    api.add_resource(AnnotationsObsAPI, "/annotations/obs")
    api.add_resource(AnnotationsVarAPI, "/annotations/var")
    api.add_resource(DataVarAPI, "/data/var")
    # Computation routes
    api.add_resource(DiffExpObsAPI, "/diffexp/obs")
    api.add_resource(LayoutObsAPI, "/layout/obs")
    return api


class ServerMulti:
    def __init__(self, app_config, annotations):

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

        resources = get_api_resources()
        self.app.register_blueprint(web_bp)
        self.app.register_blueprint(resources.blueprint)
        self.app.add_url_rule("/<dataset>/", 'index', dataset_index)
        self.app.add_url_rule("/<dataset>/static/<path:therest>", "static_redirect", static_redirect)
        self.app.add_url_rule("/<dataset>/favicon.png", "favicon_redirect", favicon_redirect)

        self.app.app_config = app_config
        self.app.annotations = annotations
