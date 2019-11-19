from http import HTTPStatus
import warnings
from uuid import uuid4
import re

from flask import Blueprint, current_app, jsonify, make_response, request, session
from flask_restful import Api, Resource
from server import __version__ as cellxgene_version
from anndata import __version__ as anndata_version

from server.app.util.constants import (
    Axis,
    DiffExpMode,
    JSON_NaN_to_num_warning_msg,
    CXGUID,
    CXG_ANNO_COLLECTION
)
from server.app.util.errors import (
    FilterError,
    InteractiveError,
    JSONEncodingValueError,
    PrepareError,
    DisabledFeatureError,
)


class SchemaAPI(Resource):
    def get(self):
        cxguid = get_userid(session)
        anno_collection = get_anno_collection(session)
        return make_response(
            jsonify({"schema": current_app.data.get_schema(uid=cxguid, collection=anno_collection)}), HTTPStatus.OK
        )


class ConfigAPI(Resource):
    def get(self):
        cxguid = get_userid(session)
        anno_collection = get_anno_collection(session)
        config = {
            "config": {
                "features": [
                    {
                        "method": "POST",
                        "path": "/cluster/",
                        **current_app.data.features["cluster"],
                    },
                    {
                        "method": "POST",
                        "path": "/layout/obs",
                        **current_app.data.features["layout"]["obs"],
                    },
                    {
                        "method": "POST",
                        "path": "/layout/var",
                        **current_app.data.features["layout"]["var"],
                    },
                    {
                        "method": "POST",
                        "path": "/diffexp/",
                        **current_app.data.features["diffexp"],
                    },
                ],
                "displayNames": {
                    "engine": f"cellxgene Scanpy engine version ",
                    "dataset": current_app.config["DATASET_TITLE"],
                },
                "links": {
                    "about-dataset": current_app.config["ABOUT_DATASET"]
                },
                "parameters": {
                    **current_app.data.get_config_parameters(uid=cxguid, collection=anno_collection)
                },
                "library_versions": {
                    "cellxgene": cellxgene_version,
                    "anndata": anndata_version
                }
            }
        }

        return make_response(jsonify(config), HTTPStatus.OK)


class AnnotationsObsAPI(Resource):
    def get(self):
        fields = request.args.getlist("annotation-name", None)
        preferred_mimetype = request.accept_mimetypes.best_match(
            ["application/octet-stream"]
        )
        cxguid = get_userid(session)
        anno_collection = get_anno_collection(session)
        try:
            if preferred_mimetype == "application/octet-stream":
                fbs = current_app.data.annotation_to_fbs_matrix("obs", fields, uid=cxguid, collection=anno_collection)
                return make_response(fbs,
                                     HTTPStatus.OK,
                                     {"Content-Type": "application/octet-stream"})
            else:
                return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
        except KeyError:
            return make_response(f"Error bad key in {fields}", HTTPStatus.BAD_REQUEST)
        except ValueError as e:
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)

    def put(self):
        cxguid = get_userid(session)
        anno_collection = request.args.get("annotation-collection-name", default=None)
        if anno_collection is not None:
            if not is_safe_collection_name(anno_collection):
                return make_response(f"Error, bad annotation collection name", HTTPStatus.BAD_REQUEST)
            set_anno_collection(session, anno_collection)
        else:
            anno_collection = get_anno_collection(session)

        try:
            fbs = request.get_data()
            res = current_app.data.annotation_put_fbs("obs", fbs, uid=cxguid, collection=anno_collection)
            return make_response(
                res, HTTPStatus.OK, {"Content-Type": "application/json"}
            )
        except (ValueError, DisabledFeatureError, KeyError) as e:
            return make_response(str(e), HTTPStatus.BAD_REQUEST)
        except Exception as e:
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


class AnnotationsVarAPI(Resource):
    def get(self):
        fields = request.args.getlist("annotation-name", None)
        preferred_mimetype = request.accept_mimetypes.best_match(
            ["application/octet-stream"]
        )
        try:
            if preferred_mimetype == "application/octet-stream":
                return make_response(current_app.data.annotation_to_fbs_matrix("var", fields),
                                     HTTPStatus.OK,
                                     {"Content-Type": "application/octet-stream"})
            else:
                return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
        except KeyError:
            return make_response(f"Error bad key in {fields}", HTTPStatus.BAD_REQUEST)
        except ValueError as e:
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


class DataVarAPI(Resource):
    def put(self):
        preferred_mimetype = request.accept_mimetypes.best_match(
            ["application/octet-stream"]
        )
        try:
            if preferred_mimetype == "application/octet-stream":
                filter_json = request.get_json()
                filter = filter_json["filter"] if filter_json else None
                return make_response(
                    current_app.data.data_frame_to_fbs_matrix(
                        filter, axis=Axis.VAR
                    ),
                    HTTPStatus.OK,
                    {"Content-Type": "application/octet-stream"})
            else:
                return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
        except FilterError as e:
            return make_response(e.message, HTTPStatus.BAD_REQUEST)
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
            return make_response(
                f"Error: invalid mode option {args['mode']}", HTTPStatus.BAD_REQUEST
            )
        # Validate filters
        if mode == DiffExpMode.VAR_FILTER or "varFilter" in args:
            # not NOT_IMPLEMENTED
            return make_response(
                "mode=varfilter not implemented", HTTPStatus.NOT_IMPLEMENTED
            )
        if mode == DiffExpMode.TOP_N and "count" not in args:
            return make_response(
                "mode=topN requires a count parameter", HTTPStatus.BAD_REQUEST
            )

        if "set1" not in args:
            return make_response("set1 is required.", HTTPStatus.BAD_REQUEST)
        if Axis.VAR in args["set1"]["filter"]:
            return make_response(
                "Var filter not allowed for set1", HTTPStatus.BAD_REQUEST
            )
        # set2
        if "set2" not in args:
            return make_response(
                "Set2 as inverse of set1 is not implemented", HTTPStatus.NOT_IMPLEMENTED
            )
        if Axis.VAR in args["set2"]["filter"]:
            return make_response(
                "Var filter not allowed for set2", HTTPStatus.BAD_REQUEST
            )

        set1_filter = args["set1"]["filter"]
        set2_filter = args.get("set2", {"filter": {}})["filter"]

        # TODO: implement varfilter mode

        # mode=topN
        count = args.get("count", None)
        try:
            diffexp = current_app.data.diffexp_topN(
                set1_filter,
                set2_filter,
                count,
                current_app.data.features["diffexp"]["interactiveLimit"],
            )
            return make_response(
                diffexp, HTTPStatus.OK, {"Content-Type": "application/json"}
            )
        except (ValueError, FilterError) as e:
            return make_response(e.message, HTTPStatus.BAD_REQUEST)
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
        preferred_mimetype = request.accept_mimetypes.best_match(
            ["application/octet-stream"]
        )
        try:
            if preferred_mimetype == "application/octet-stream":
                return make_response(current_app.data.layout_to_fbs_matrix(),
                                     HTTPStatus.OK,
                                     {"Content-Type": "application/octet-stream"})
            else:
                return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
        except PrepareError as e:
            return make_response(e.message, HTTPStatus.INTERNAL_SERVER_ERROR)
        except ValueError as e:
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


def get_userid(ss):
    if CXGUID not in ss:
        ss[CXGUID] = uuid4().hex
        ss.permanent = True
    return ss[CXGUID]


def get_anno_collection(ss):
    collection = ss[CXG_ANNO_COLLECTION] if CXG_ANNO_COLLECTION in ss else None
    return collection


def set_anno_collection(ss, name):
    ss[CXG_ANNO_COLLECTION] = name
    ss.permanent = True


def is_safe_collection_name(name):
    """
    return true if this is a safe collection name

    this is ultra convervative. If we want to allow full legal file name syntax,
    we could look at modules like `pathvalidate`
    """
    if name is None:
        return False
    return re.match(r'^\w+$', name) is not None


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
