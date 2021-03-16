import datetime
import logging
from functools import wraps
from http import HTTPStatus

from flask import (
    Flask,
    current_app,
    make_response,
    render_template,
    Blueprint,
    request,
    send_from_directory,
)
from flask_restful import Api, Resource

import backend.server.common.rest as common_rest
from backend.common_utils.errors import DatasetAccessError, RequestException
from backend.server.common.health import health_check
from backend.common_utils.utils import Float32JSONEncoder

webbp = Blueprint("webapp", "backend.server.common.web", template_folder="templates")


@webbp.route("/", methods=["GET"])
def dataset_index():
    app_config = current_app.app_config

    dataset_config = app_config.get_dataset_config()
    scripts = dataset_config.app__scripts
    inline_scripts = dataset_config.app__inline_scripts

    try:
        args = {"SCRIPTS": scripts, "INLINE_SCRIPTS": inline_scripts}
        return render_template("index.html", **args)

    except DatasetAccessError as e:
        return common_rest.abort_and_log(
            e.status_code, f"Invalid dataset: {e.message}", loglevel=logging.INFO, include_exc_info=True
        )


@webbp.errorhandler(RequestException)
def handle_request_exception(error):
    return common_rest.abort_and_log(error.status_code, error.message, loglevel=logging.INFO, include_exc_info=True)


def requires_authentication(func):
    @wraps(func)
    def wrapped_function(self, *args, **kwargs):
        auth = current_app.auth
        if auth.is_user_authenticated():
            return func(self, *args, **kwargs)
        else:
            return make_response("not authenticated", HTTPStatus.UNAUTHORIZED)

    return wrapped_function


def rest_get_data_adaptor(func):
    @wraps(func)
    def wrapped_function(self):
        try:
            return func(self, current_app.data_adaptor)
        except DatasetAccessError as e:
            return common_rest.abort_and_log(
                e.status_code, f"Invalid dataset: {e.message}", loglevel=logging.INFO, include_exc_info=True
            )

    return wrapped_function


class HealthAPI(Resource):
    def get(self):
        config = current_app.app_config
        return health_check(config)


class SchemaAPI(Resource):
    # TODO @mdunitz separate dataset schema and user schema
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.schema_get(data_adaptor)


class ConfigAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.config_get(current_app.app_config, data_adaptor)


class UserInfoAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.userinfo_get(current_app.app_config, data_adaptor)


class AnnotationsObsAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.annotations_obs_get(request, data_adaptor)

    @requires_authentication
    @rest_get_data_adaptor
    def put(self, data_adaptor):
        return common_rest.annotations_obs_put(request, data_adaptor)


class AnnotationsVarAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.annotations_var_get(request, data_adaptor)


class DataVarAPI(Resource):
    @rest_get_data_adaptor
    def put(self, data_adaptor):
        return common_rest.data_var_put(request, data_adaptor)

    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.data_var_get(request, data_adaptor)


class ColorsAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.colors_get(data_adaptor)


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


class GenesetsAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.genesets_get(request, data_adaptor)

    @requires_authentication
    @rest_get_data_adaptor
    def put(self, data_adaptor):
        return common_rest.genesets_put(request, data_adaptor)


class GenesetSummaryAPI(Resource):
    @rest_get_data_adaptor
    def get(self, data_adaptor):
        return common_rest.geneset_summary_get(request, data_adaptor)


def get_api_base_resources(bp_base):
    """Add resources that are accessed from the api url"""
    api = Api(bp_base)

    # Diagnostics routes
    api.add_resource(HealthAPI, "/health")
    return api


def get_api_dataroot_resources(bp_dataroot):
    """Add resources that refer to a dataset"""
    api = Api(bp_dataroot)

    def add_resource(resource, url):
        """convenience function to make the outer function less verbose"""
        api.add_resource(resource, url)

    # Initialization routes
    add_resource(SchemaAPI, "/schema")
    add_resource(ConfigAPI, "/config")
    add_resource(UserInfoAPI, "/userinfo")
    # Data routes
    add_resource(AnnotationsObsAPI, "/annotations/obs")
    add_resource(AnnotationsVarAPI, "/annotations/var")
    add_resource(DataVarAPI, "/data/var")
    add_resource(GenesetsAPI, "/genesets")
    add_resource(GenesetSummaryAPI, "/geneset_summary")
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
        self.app = Flask(__name__, static_folder=None)
        self._before_adding_routes(self.app, app_config)
        self.app.json_encoder = Float32JSONEncoder
        server_config = app_config.server_config

        # enable session data
        self.app.permanent_session_lifetime = datetime.timedelta(days=50 * 365)

        # Config
        secret_key = server_config.app__flask_secret_key
        self.app.config.update(SECRET_KEY=secret_key)

        self.app.register_blueprint(webbp)

        api_version = "/api/v0.2"
        api_path = "/"

        bp_base = Blueprint("bp_base", __name__, url_prefix=api_path)
        base_resources = get_api_base_resources(bp_base)
        self.app.register_blueprint(base_resources.blueprint)

        bp_api = Blueprint("api", __name__, url_prefix=f"{api_path}{api_version}")
        resources = get_api_dataroot_resources(bp_api)
        self.app.register_blueprint(resources.blueprint)
        self.app.add_url_rule(
            "/static/<path:filename>",
            "static_assets",
            view_func=lambda filename: send_from_directory("../common/web/static", filename),
            methods=["GET"],
        )

        self.app.data_adaptor = server_config.data_adaptor
        self.app.app_config = app_config

        auth = server_config.auth
        self.app.auth = auth
        if auth.requires_client_login():
            auth.add_url_rules(self.app)
        auth.complete_setup(self.app)
