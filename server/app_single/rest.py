from flask import Blueprint, current_app, request
from flask_restful import Api, Resource
import server.common.rest as common_rest


class SchemaAPI(Resource):
    def get(self):
        return common_rest.schema_get(current_app.data, current_app.annotations)


class ConfigAPI(Resource):
    def get(self):
        return common_rest.config_get(
            current_app.app_config, current_app.data, current_app.annotations)


class AnnotationsObsAPI(Resource):
    def get(self):
        return common_rest.annotations_obs_get(request, current_app.data, current_app.annotations)

    def put(self):

        return common_rest.annotations_obs_put(
            request, current_app.data, current_app.annotations)


class AnnotationsVarAPI(Resource):
    def get(self):
        return common_rest.annotations_var_get(request, current_app.data, current_app.annotations)


class DataVarAPI(Resource):
    def put(self):
        return common_rest.data_var_put(request, current_app.data)


class DiffExpObsAPI(Resource):
    def post(self):
        return common_rest.diffexp_obs_post(request, current_app.data)


class LayoutObsAPI(Resource):
    def get(self):
        return common_rest.layout_obs_get(request, current_app.data)


def get_api_resources():
    bp_api = Blueprint("api", __name__, url_prefix="/api/v0.2")
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
