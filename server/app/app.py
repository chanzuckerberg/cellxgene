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
from server.common.errors import DatasetAccessError
from server.common.utils import path_join, Float32JSONEncoder
from server.common.data_locator import DataLocator
from server.data_common.matrix_loader import MatrixDataLoader, MatrixDataType

from functools import wraps

webbp = Blueprint("webapp", "server.common.web", template_folder="templates")


@webbp.route("/")
def dataset_index(dataset=None):
    config = current_app.app_config
    if dataset is None:
        if config.datapath:
            location = config.datapath
        else:
            return dataroot_index()
    else:
        location = path_join(config.dataroot, dataset)

    scripts = config.scripts

    try:
        cache_manager = current_app.matrix_data_cache_manager
        with cache_manager.data_adaptor(location, config) as data_adaptor:
            dataset_title = config.get_title(data_adaptor)
            return render_template("index.html", datasetTitle=dataset_title, SCRIPTS=scripts)
    except DatasetAccessError as e:
        return make_response(f"Invalid dataset {dataset}: {str(e)}", HTTPStatus.BAD_REQUEST)


@webbp.route("/favicon.png")
def favicon():
    return send_from_directory(os.path.join(webbp.root_path, "static/img/"), "favicon.png")


def get_data_adaptor(dataset=None):
    config = current_app.app_config

    if dataset is None:
        datapath = config.datapath
    else:
        datapath = path_join(config.dataroot, dataset)
        # path_join returns a normalized path.  Therefore it is
        # sufficient to check that the datapath starts with the
        # dataroot to determine that the datapath is under the dataroot.
        if not datapath.startswith(config.dataroot):
            raise DatasetAccessError("Invalid dataset {dataset}")

    if datapath is None:
        return make_response("Dataset must be supplied", HTTPStatus.BAD_REQUEST)

    cache_manager = current_app.matrix_data_cache_manager
    return cache_manager.data_adaptor(datapath, config)


def rest_get_data_adaptor(func):
    @wraps(func)
    def wrapped_function(self, dataset=None):
        try:
            with get_data_adaptor(dataset) as data_adaptor:
                return func(self, data_adaptor)
        except DatasetAccessError as e:
            return make_response(f"Invalid dataset {dataset}: {str(e)}", HTTPStatus.BAD_REQUEST)

    return wrapped_function


def static_redirect(dataset, therest):
    """ redirect all static requests to the standard location """
    return redirect(f"/static/{therest}", code=301)


def favicon_redirect(dataset):
    """ redirect favicon to static dir """
    return redirect("/static/favicon.png", code=301)


def dataroot_index():
    # FIXME with a splash screen that includes a listing of all the datasets.
    # or perhaps a login screen if this is a hosted environment
    data = "<H1>Welcome to cellxgene</H1>"

    # the following is just for demo purposes...
    try:
        config = current_app.app_config
        locator = DataLocator(config.dataroot)
        datasets = []
        for fname in locator.ls():
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
    except Exception:
        pass

    return make_response(data)


class SchemaAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.schema_get(data_adaptor, current_app.annotations)


class ConfigAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.config_get(current_app.app_config, data_adaptor, current_app.annotations)


class AnnotationsObsAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.annotations_obs_get(request, data_adaptor, current_app.annotations)

    @rest_get_data_adaptor
    def put(self, data_adaptor):
        return common_rest.annotations_obs_put(request, data_adaptor, current_app.annotations)


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

    @rest_get_data_adaptor
    def put(self, data_adaptor):
        return common_rest.layout_obs_put(request, data_adaptor)


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
    def __init__(self, matrix_data_cache_manager, annotations, app_config):

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
        if app_config.datapath:
            bp_api = Blueprint("api", __name__, url_prefix=api_version)
            resources = get_api_resources(bp_api)
            self.app.register_blueprint(resources.blueprint)

        else:
            # NOTE:  These routes only allow the dataset to be in the directory
            # of the dataroot, and not a subdirectory.  We may want to change
            # the route format at some point
            bp_api = Blueprint("api_dataset", __name__, url_prefix="/<dataset>" + api_version)
            resources = get_api_resources(bp_api)
            self.app.register_blueprint(resources.blueprint)
            self.app.add_url_rule("/<dataset>/", "dataset_index", dataset_index)
            self.app.add_url_rule("/<dataset>/static/<path:therest>", "static_redirect", static_redirect)
            self.app.add_url_rule("/<dataset>/favicon.png", "favicon_redirect", favicon_redirect)

        self.app.matrix_data_cache_manager = matrix_data_cache_manager
        self.app.annotations = annotations
        self.app.app_config = app_config
