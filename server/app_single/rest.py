from http import HTTPStatus
import warnings

from flask import Blueprint, current_app, jsonify, make_response, request
from flask_restful import Api, Resource

from server.common.constants import Axis, DiffExpMode, JSON_NaN_to_num_warning_msg
from server.common.errors import (
    FilterError,
    InteractiveError,
    JSONEncodingValueError,
    PrepareError,
    DisabledFeatureError,
)

import json
from server.common.utils import get_schema, annotation_put_fbs


class SchemaAPI(Resource):
    def get(self):
        schema = get_schema(current_app.data, current_app.annotations)
        return make_response(
            jsonify({"schema": schema}), HTTPStatus.OK
        )


class ConfigAPI(Resource):
    def get(self):
        config = current_app.app_config.get_config(
            current_app.data,
            current_app.annotations)
        return make_response(make_response(jsonify(config), HTTPStatus.OK))


class AnnotationsObsAPI(Resource):
    def get(self):
        fields = request.args.getlist("annotation-name", None)
        preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
        try:
            if preferred_mimetype == "application/octet-stream":
                labels = None
                if current_app.annotations:
                    labels = current_app.annotations.read_labels()
                fbs = current_app.data.annotation_to_fbs_matrix(Axis.OBS, fields, labels)
                return make_response(fbs, HTTPStatus.OK, {"Content-Type": "application/octet-stream"})
            else:
                return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
        except KeyError:
            return make_response(f"Error bad key in {fields}", HTTPStatus.BAD_REQUEST)
        except ValueError as e:
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)
        except Exception as e:
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)

    def put(self):
        anno_collection = request.args.get("annotation-collection-name", default=None)
        annotations = current_app.annotations
        if annotations is None:
            return make_response("Error, annotations are not configured", HTTPStatus.BAD_REQUEST)

        if anno_collection is not None:
            if not annotations.is_safe_collection_name(anno_collection):
                return make_response(f"Error, bad annotation collection name", HTTPStatus.BAD_REQUEST)
            annotations.set_collection(anno_collection)

        try:
            fbs = request.get_data()
            annotation_put_fbs(fbs, current_app.data, current_app.annotations)
            res = json.dumps({"status": "OK"})
            return make_response(res, HTTPStatus.OK, {"Content-Type": "application/json"})
        except (ValueError, DisabledFeatureError, KeyError) as e:
            return make_response(str(e), HTTPStatus.BAD_REQUEST)
        except Exception as e:
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


class AnnotationsVarAPI(Resource):
    def get(self):
        fields = request.args.getlist("annotation-name", None)

        preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
        try:
            if preferred_mimetype == "application/octet-stream":
                labels = None
                if current_app.annotations is not None:
                    labels = current_app.annotations.read_labels()
                return make_response(
                    current_app.data.annotation_to_fbs_matrix(Axis.VAR, fields, labels),
                    HTTPStatus.OK,
                    {"Content-Type": "application/octet-stream"},
                )
            else:
                return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
        except KeyError:
            return make_response(f"Error bad key in {fields}", HTTPStatus.BAD_REQUEST)
        except ValueError as e:
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)
        except Exception as e:
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


class DataVarAPI(Resource):
    def put(self):
        preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
        try:
            if preferred_mimetype == "application/octet-stream":
                filter_json = request.get_json()
                filter = filter_json["filter"] if filter_json else None
                return make_response(
                    current_app.data.data_frame_to_fbs_matrix(filter, axis=Axis.VAR),
                    HTTPStatus.OK,
                    {"Content-Type": "application/octet-stream"},
                )
            else:
                return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
        except FilterError as e:
            return make_response(str(e), HTTPStatus.BAD_REQUEST)
        except ValueError as e:
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


class DiffExpObsAPI(Resource):
    def post(self):
        args = request.get_json()
        # confirm mode is present and legal
        try:
            mode = DiffExpMode(args["mode"])
        except KeyError:
            return make_response("Error: mode is required", HTTPStatus.BAD_REQUEST)
        except ValueError:
            return make_response(f"Error: invalid mode option {args['mode']}", HTTPStatus.BAD_REQUEST)
        # Validate filters
        if mode == DiffExpMode.VAR_FILTER or "varFilter" in args:
            # not NOT_IMPLEMENTED
            return make_response("mode=varfilter not implemented", HTTPStatus.NOT_IMPLEMENTED)
        if mode == DiffExpMode.TOP_N and "count" not in args:
            return make_response("mode=topN requires a count parameter", HTTPStatus.BAD_REQUEST)

        if "set1" not in args:
            return make_response("set1 is required.", HTTPStatus.BAD_REQUEST)
        if Axis.VAR in args["set1"]["filter"]:
            return make_response("Var filter not allowed for set1", HTTPStatus.BAD_REQUEST)
        # set2
        if "set2" not in args:
            return make_response("Set2 as inverse of set1 is not implemented", HTTPStatus.NOT_IMPLEMENTED)
        if Axis.VAR in args["set2"]["filter"]:
            return make_response("Var filter not allowed for set2", HTTPStatus.BAD_REQUEST)

        set1_filter = args["set1"]["filter"]
        set2_filter = args.get("set2", {"filter": {}})["filter"]

        # TODO: implement varfilter mode

        # mode=topN
        count = args.get("count", None)
        try:
            interactive_limit = current_app.data.get_features()["diffexp"].interactiveLimit
            diffexp = current_app.data.diffexp_topN(set1_filter, set2_filter, count, interactive_limit)
            return make_response(diffexp, HTTPStatus.OK, {"Content-Type": "application/json"})
        except (ValueError, FilterError) as e:
            return make_response(str(e), HTTPStatus.BAD_REQUEST)
        except InteractiveError:
            return make_response("Non-interactive request", HTTPStatus.FORBIDDEN)
        except JSONEncodingValueError as e:
            # JSON encoding failure, usually due to bad data
            warnings.warn(JSON_NaN_to_num_warning_msg)
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)
        except ValueError as e:
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


class LayoutObsAPI(Resource):
    def get(self):
        preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
        try:
            if preferred_mimetype == "application/octet-stream":
                return make_response(
                    current_app.data.layout_to_fbs_matrix(), HTTPStatus.OK, {"Content-Type": "application/octet-stream"}
                )
            else:
                return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
        except PrepareError as e:
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)
        except ValueError as e:
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


def get_api_resources():
    bp = Blueprint("api", __name__, url_prefix="/api/v0.2")
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
