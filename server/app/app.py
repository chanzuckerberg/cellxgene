import datetime
import logging

from flask import Flask, redirect, current_app, make_response, render_template, abort
from flask import Blueprint, request
from flask_restful import Api, Resource
from server_timing import Timing as ServerTiming

from http import HTTPStatus

import server.common.rest as common_rest
from server.common.errors import DatasetAccessError
from server.common.utils import path_join, Float32JSONEncoder
from server.common.data_locator import DataLocator
from server.common.health import health_check
from server.data_common.matrix_loader import MatrixDataLoader

from functools import wraps

webbp = Blueprint("webapp", "server.common.web", template_folder="templates")

ONE_WEEK = 7 * 24 * 60 * 60


def _cache_control(always, **cache_kwargs):
    """
    Used to easily manage cache control headers on responses.
    See Werkzeug for attributes that can be set, eg, no_cache, private, max_age, etc.
    https://werkzeug.palletsprojects.com/en/1.0.x/datastructures/#werkzeug.datastructures.ResponseCacheControl
    """

    def inner_cache_control(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            response = make_response(f(*args, **kwargs))
            if not always and not current_app.app_config.server__generate_cache_control_headers:
                return response
            if response.status_code >= 400:
                return response
            for k, v in cache_kwargs.items():
                setattr(response.cache_control, k, v)
            return response

        return wrapper

    return inner_cache_control


def cache_control(**cache_kwargs):
    """ configu driven """
    return _cache_control(False, **cache_kwargs)


def cache_control_always(**cache_kwargs):
    """ always generate headers, regardless of the config """
    return _cache_control(True, **cache_kwargs)


@webbp.route("/", methods=["GET"])
@cache_control(public=True, max_age=ONE_WEEK)
def dataset_index(url_dataroot=None, dataset=None):
    config = current_app.app_config
    if dataset is None:
        if config.single_dataset__datapath:
            location = config.single_dataset__datapath
        else:
            return dataroot_index()
    else:
        dataroot = config.multi_dataset__dataroot.get(url_dataroot)
        if dataroot is None:
            abort(HTTPStatus.NOT_FOUND)
        location = path_join(dataroot, dataset)

    scripts = config.server__scripts
    inline_scripts = config.server__inline_scripts

    try:
        cache_manager = current_app.matrix_data_cache_manager
        with cache_manager.data_adaptor(location, config) as data_adaptor:
            dataset_title = config.get_title(data_adaptor)
            return render_template(
                "index.html", datasetTitle=dataset_title, SCRIPTS=scripts, INLINE_SCRIPTS=inline_scripts
            )
    except DatasetAccessError:
        return common_rest.abort_and_log(
            HTTPStatus.BAD_REQUEST, f"Invalid dataset {dataset}", loglevel=logging.INFO, include_exc_info=True
        )


@webbp.route("/health", methods=["GET"])
@cache_control_always(no_store=True)
def health():
    config = current_app.app_config
    return health_check(config)


def get_data_adaptor(url_dataroot=None, dataset=None):
    config = current_app.app_config

    if dataset is None:
        datapath = config.single_dataset__datapath
    else:
        dataroot = config.multi_dataset__dataroot.get(url_dataroot)
        if dataroot is None:
            raise DatasetAccessError(f"Invalid dataset {url_dataroot}/{dataset}")
        datapath = path_join(dataroot, dataset)
        # path_join returns a normalized path.  Therefore it is
        # sufficient to check that the datapath starts with the
        # dataroot to determine that the datapath is under the dataroot.
        if not datapath.startswith(dataroot):
            raise DatasetAccessError("Invalid dataset {url_dataroot}/{dataset}")

    if datapath is None:
        return common_rest.abort_and_log(HTTPStatus.BAD_REQUEST, "Invalid dataset NONE", loglevel=logging.INFO)

    cache_manager = current_app.matrix_data_cache_manager
    return cache_manager.data_adaptor(datapath, config)


def rest_get_data_adaptor(func):
    @wraps(func)
    def wrapped_function(self, dataset=None):
        try:
            with get_data_adaptor(self.url_dataroot, dataset) as data_adaptor:
                return func(self, data_adaptor)
        except DatasetAccessError:
            return common_rest.abort_and_log(
                HTTPStatus.BAD_REQUEST, f"Invalid dataset {dataset}", loglevel=logging.INFO, include_exc_info=True
            )

    return wrapped_function


def dataroot_test_index():
    # the following index page is meant for testing/debugging purposes
    data = '<!doctype html><html lang="en">'
    data += "<head><title>Hosted Cellxgene</title></head>"
    data += "<body><H1>Welcome to cellxgene</H1>"

    config = current_app.app_config
    datasets = []
    for url_dataroot, dataroot in config.multi_dataset__dataroot.items():
        locator = DataLocator(dataroot, region_name=config.data_locator__s3__region_name)
        for fname in locator.ls():
            location = path_join(dataroot, fname)
            try:
                MatrixDataLoader(location, app_config=config)
                datasets.append((url_dataroot, fname))
            except DatasetAccessError:
                # skip over invalid datasets
                pass

    data += "<br/>Select one of these datasets...<br/>"
    data += "<ul>"
    datasets.sort()
    for url_dataroot, dataset in datasets:
        data += f"<li><a href={url_dataroot}/{dataset}>{dataset}</a></li>"
    data += "</ul>"
    data += "</body></html>"

    return make_response(data)


def dataroot_index():
    # Handle the base url for the cellxgene server when running in multi dataset mode
    config = current_app.app_config
    if not config.multi_dataset__index:
        abort(HTTPStatus.NOT_FOUND)
    elif config.multi_dataset__index is True:
        return dataroot_test_index()
    else:
        return redirect(config.multi_dataset__index)


class DatasetResource(Resource):
    """Base class for all Resources that act on datasets."""

    def __init__(self, url_dataroot):
        super().__init__()
        self.url_dataroot = url_dataroot


class SchemaAPI(DatasetResource):
    @cache_control(public=True, max_age=ONE_WEEK)
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.schema_get(data_adaptor, current_app.annotations)


class ConfigAPI(DatasetResource):
    @cache_control(public=True, max_age=ONE_WEEK)
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.config_get(current_app.app_config, data_adaptor, current_app.annotations)


class AnnotationsObsAPI(DatasetResource):
    @cache_control(public=True, max_age=ONE_WEEK)
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.annotations_obs_get(request, data_adaptor, current_app.annotations)

    @cache_control(no_store=True)
    @rest_get_data_adaptor
    def put(self, data_adaptor):
        return common_rest.annotations_obs_put(request, data_adaptor, current_app.annotations)


class AnnotationsVarAPI(DatasetResource):
    @cache_control(public=True, max_age=ONE_WEEK)
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.annotations_var_get(request, data_adaptor, current_app.annotations)


class DataVarAPI(DatasetResource):
    @cache_control(no_store=True)
    @rest_get_data_adaptor
    def put(self, data_adaptor):
        return common_rest.data_var_put(request, data_adaptor)

    @cache_control(public=True, max_age=ONE_WEEK)
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.data_var_get(request, data_adaptor)


class ColorsAPI(DatasetResource):
    @cache_control(public=True, max_age=ONE_WEEK)
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.colors_get(data_adaptor)


class DiffExpObsAPI(DatasetResource):
    @cache_control(no_store=True)
    @rest_get_data_adaptor
    def post(self, data_adaptor):
        return common_rest.diffexp_obs_post(request, data_adaptor)


class LayoutObsAPI(DatasetResource):
    @cache_control(public=True, max_age=ONE_WEEK)
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.layout_obs_get(request, data_adaptor)

    @cache_control(no_store=True)
    @rest_get_data_adaptor
    def put(self, data_adaptor):
        return common_rest.layout_obs_put(request, data_adaptor)


def get_api_resources(bp_api, url_dataroot=None):
    api = Api(bp_api)

    def add_resource(resource, url):
        """convenience function to make the outer function less verbose"""
        api.add_resource(resource, url, resource_class_args=(url_dataroot,))

    # Initialization routes
    add_resource(SchemaAPI, "/schema")
    add_resource(ConfigAPI, "/config")
    # Data routes
    add_resource(AnnotationsObsAPI, "/annotations/obs")
    add_resource(AnnotationsVarAPI, "/annotations/var")
    add_resource(DataVarAPI, "/data/var")
    # Display routes
    add_resource(ColorsAPI, "/colors")
    # Computation routes
    add_resource(DiffExpObsAPI, "/diffexp/obs")
    add_resource(LayoutObsAPI, "/layout/obs")
    return api


class Server:
    @staticmethod
    def _before_adding_routes(app, app_config):
        """ will be called before routes are added, during __init__.  Subclass protocol """
        pass

    def __init__(self, app_config):
        self.app = Flask(__name__, static_folder="../common/web/static")
        self._before_adding_routes(self.app, app_config)
        self.app.json_encoder = Float32JSONEncoder
        if app_config.server__server_timing_headers:
            ServerTiming(self.app, force_debug=True)

        # enable session data
        self.app.permanent_session_lifetime = datetime.timedelta(days=50 * 365)

        # Config
        secret_key = app_config.server__flask_secret_key
        self.app.config.update(SECRET_KEY=secret_key)

        self.app.register_blueprint(webbp)

        api_version = "/api/v0.2"
        if app_config.single_dataset__datapath:
            bp_api = Blueprint("api", __name__, url_prefix=api_version)
            resources = get_api_resources(bp_api)
            self.app.register_blueprint(resources.blueprint)

        else:
            # NOTE:  These routes only allow the dataset to be in the directory
            # of the dataroot, and not a subdirectory.  We may want to change
            # the route format at some point
            for url_dataroot in app_config.multi_dataset__dataroot.keys():
                bp_api = Blueprint(
                    f"api_dataset_{url_dataroot}", __name__, url_prefix=f"/{url_dataroot}/<dataset>" + api_version
                )
                resources = get_api_resources(bp_api, url_dataroot)
                self.app.register_blueprint(resources.blueprint)
                self.app.add_url_rule(
                    f"/{url_dataroot}/<dataset>/",
                    f"dataset_index_{url_dataroot}",
                    lambda dataset: dataset_index(url_dataroot, dataset),
                    methods=["GET"],
                )
        self.app.matrix_data_cache_manager = app_config.matrix_data_cache_manager
        self.app.annotations = app_config.user_annotations
        self.app.app_config = app_config
