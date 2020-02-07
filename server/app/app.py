import os
import datetime

from flask import Flask, redirect, current_app, make_response, render_template
from flask import Blueprint, request, send_from_directory
from flask_caching import Cache
from flask_compress import Compress
from flask_cors import CORS
from flask_restful import Api, Resource

from http import HTTPStatus

import server.common.rest as common_rest
from server.common.utils import path_join, Float32JSONEncoder
from server.data_common.matrix_loader import MatrixDataLoader, MatrixDataType, MatrixDataCacheManager

from functools import wraps

webbp = Blueprint("webapp", "server.common.web", template_folder="templates")


@webbp.route("/")
def index():
    if current_app.data is None:
        return dataroot_index()

    dataset_title = current_app.app_config.title
    scripts = current_app.app_config.scripts
    return render_template("index.html", datasetTitle=dataset_title, SCRIPTS=scripts)


@webbp.route("/favicon.png")
def favicon():
    return send_from_directory(os.path.join(webbp.root_path, "static/img/"), "favicon.png")


def get_data_adaptor(dataset):
    config = current_app.app_config
    location = path_join(config.dataroot, dataset)
    matrix_data_loader = MatrixDataLoader(location)
    matrix_data_loader.pre_load_validation()
    data = matrix_data_loader.open(current_app.app_config)
    return data


def rest_get_data_adaptor(func):
    @wraps(func)
    def wrapped_function(self, dataset=None):
        try:
            if dataset is None:
                # use the default dataset
                data_adaptor = current_app.data
                if data_adaptor is None:
                    return make_response("Dataset must be supplied", HTTPStatus.BAD_REQUEST)
                return func(self, data_adaptor)

            config = current_app.app_config
            location = path_join(config.dataroot, dataset)
            cache_manager = current_app.matrix_data_cache_manager
            with cache_manager.data_adaptor(location, config) as data_adaptor:
                return func(self, data_adaptor)
        except RuntimeError as e:
            return make_response(f"Invalid dataset {dataset}: {str(e)}", HTTPStatus.BAD_REQUEST)
    return wrapped_function


def dataset_index(dataset):
    config = current_app.app_config
    location = path_join(config.dataroot, dataset)
    matrix_data_loader = MatrixDataLoader(location)
    try:
        matrix_data_loader.pre_load_validation()
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


def dataroot_index():
    # FIXME with a splash screen that includes a listing of all the datasets.
    # or perhaps a login screen if this is a hosted environment
    data = "<H1>Welcome to cellxgene</H1>"

    # the following is just for demo purposes...
    config = current_app.app_config
    datasets = []
    for fname in os.listdir(config.dataroot):
        location = path_join(config.dataroot, fname)
        matrix_data_loader = MatrixDataLoader(location)
        if matrix_data_loader.etype != MatrixDataType.UNKNOWN:
            datasets.append(fname)

    data += "<br/>Select one of these datasets...<br/>"
    data += "<ul>"
    datasets.sort()
    for dataset in datasets:
        data += f"<li><a href={dataset}>{dataset}</a></li>"
    data += "</ul>"

    return make_response(data)


class SchemaAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.schema_get(data_adaptor, current_app.annotations)


class ConfigAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.config_get(
            current_app.app_config, data_adaptor, current_app.annotations)


class AnnotationsObsAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.annotations_obs_get(
            request, data_adaptor, current_app.annotations)

    @rest_get_data_adaptor
    def put(self, data_adaptor):
        return common_rest.annotations_obs_put(
            request, data_adaptor, current_app.annotations)


class AnnotationsVarAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.annotations_var_get(request, data_adaptor, current_app.annotations)


class DataVarAPI(Resource):
    @rest_get_data_adaptor
    def put(self, data_adaptor):
        return common_rest.data_var_put(request, data_adaptor)


class DiffExpObsAPI(Resource):
    @rest_get_data_adaptor
    def post(self, data_adaptor):
        return common_rest.diffexp_obs_post(request, data_adaptor)


class LayoutObsAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.layout_obs_get(request, data_adaptor)


def get_api_resources(bp_api):
    api = Api(bp_api)
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


class Server:
    def __init__(self, data, annotations, app_config):

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

        self.app.register_blueprint(webbp)

        api_version = "/api/v0.2"
        if data:
            bp_api = Blueprint("api", __name__, url_prefix=api_version)
            resources = get_api_resources(bp_api)
            self.app.register_blueprint(resources.blueprint)
            self.app.add_url_rule("/", "index", index)

        if app_config.dataroot:
            bp_api = Blueprint("api_dataset", __name__, url_prefix="/<dataset>" + api_version)
            resources = get_api_resources(bp_api)
            self.app.register_blueprint(resources.blueprint)
            self.app.add_url_rule("/<dataset>/", 'dataset_index', dataset_index)
            self.app.add_url_rule("/<dataset>/static/<path:therest>", "static_redirect", static_redirect)
            self.app.add_url_rule("/<dataset>/favicon.png", "favicon_redirect", favicon_redirect)

        self.app.data = data
        self.app.annotations = annotations
        self.app.app_config = app_config
        self.app.matrix_data_cache_manager = MatrixDataCacheManager()
