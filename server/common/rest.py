import sys
from http import HTTPStatus
import copy
import logging
from flask import make_response, jsonify, current_app, abort
from server.common.constants import Axis, DiffExpMode, JSON_NaN_to_num_warning_msg
from server.common.errors import (
    FilterError,
    JSONEncodingValueError,
    PrepareError,
    DisabledFeatureError,
)

import json
from server.data_common.fbs.matrix import decode_matrix_fbs


def abort_and_log(code, logmsg, loglevel=logging.DEBUG, include_exc_info=False):
    """
    Log the message, then abort with HTTP code. If include_exc_info is true,
    also include current exception via sys.exc_info().
    """
    if include_exc_info:
        exc_info = sys.exc_info()
    else:
        exc_info = False
    current_app.logger.log(loglevel, logmsg, exc_info=exc_info)
    # Do NOT send log message to HTTP response.
    return abort(code)


def schema_get_helper(data_adaptor, annotations):
    """helper function to gather the schema from the data source and annotations"""
    schema = data_adaptor.get_schema()
    schema = copy.deepcopy(schema)

    # add label obs annotations as needed
    if annotations is not None:
        label_schema = annotations.get_schema(data_adaptor)
        schema["annotations"]["obs"]["columns"].extend(label_schema)

    return schema


def schema_get(data_adaptor, annotations):
    schema = schema_get_helper(data_adaptor, annotations)
    return make_response(jsonify({"schema": schema}), HTTPStatus.OK)


def config_get(app_config, data_adaptor, annotations):
    config = app_config.get_config(data_adaptor, annotations)
    return make_response(jsonify(config), HTTPStatus.OK)


def annotations_obs_get(request, data_adaptor, annotations):
    fields = request.args.getlist("annotation-name", None)
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)
    try:
        labels = None
        if annotations:
            labels = annotations.read_labels(data_adaptor)
        fbs = data_adaptor.annotation_to_fbs_matrix(Axis.OBS, fields, labels)
        return make_response(fbs, HTTPStatus.OK, {"Content-Type": "application/octet-stream"})
    except KeyError as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def annotations_put_fbs_helper(data_adaptor, annotations, fbs):
    """helper function to write annotations from fbs"""
    if annotations is None:
        raise DisabledFeatureError("Writable annotations are not enabled")

    new_label_df = decode_matrix_fbs(fbs)
    if not new_label_df.empty:
        data_adaptor.check_new_labels(new_label_df)
    annotations.write_labels(new_label_df, data_adaptor)


def annotations_obs_put(request, data_adaptor, annotations):
    if annotations is None:
        return abort(HTTPStatus.NOT_IMPLEMENTED)

    anno_collection = request.args.get("annotation-collection-name", default=None)
    fbs = request.get_data()

    if anno_collection is not None:
        if not annotations.is_safe_collection_name(anno_collection):
            return abort(HTTPStatus.BAD_REQUEST, "Bad annotation collection name")
        annotations.set_collection(anno_collection)

    try:
        annotations_put_fbs_helper(data_adaptor, annotations, fbs)
        res = json.dumps({"status": "OK"})
        return make_response(res, HTTPStatus.OK, {"Content-Type": "application/json"})
    except (ValueError, DisabledFeatureError, KeyError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def annotations_var_get(request, data_adaptor, annotations):
    fields = request.args.getlist("annotation-name", None)
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)
    try:
        labels = None
        if annotations is not None:
            labels = annotations.read_labels(data_adaptor)
        return make_response(
            data_adaptor.annotation_to_fbs_matrix(Axis.VAR, fields, labels),
            HTTPStatus.OK,
            {"Content-Type": "application/octet-stream"},
        )
    except KeyError as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def data_var_put(request, data_adaptor):
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    filter_json = request.get_json()
    filter = filter_json["filter"] if filter_json else None
    try:
        return make_response(
            data_adaptor.data_frame_to_fbs_matrix(filter, axis=Axis.VAR),
            HTTPStatus.OK,
            {"Content-Type": "application/octet-stream"},
        )
    except FilterError as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def diffexp_obs_post(request, data_adaptor):
    if data_adaptor.config.disable_diffexp:
        return abort(HTTPStatus.NOT_IMPLEMENTED)

    args = request.get_json()
    try:
        # TODO: implement varfilter mode
        mode = DiffExpMode(args["mode"])

        if mode == DiffExpMode.VAR_FILTER or "varFilter" in args:
            return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, "varFilter not enabled")

        set1_filter = args.get("set1", {"filter": {}})["filter"]
        set2_filter = args.get("set2", {"filter": {}})["filter"]
        count = args.get("count", None)

        if set1_filter is None or set2_filter is None or count is None:
            return abort_and_log(HTTPStatus.BAD_REQUEST, "missing required parameter")
        if Axis.VAR in set1_filter or Axis.VAR in set2_filter:
            return abort_and_log(HTTPStatus.BAD_REQUEST, "var axis filter not enabled")

    except (KeyError, TypeError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)

    try:
        diffexp = data_adaptor.diffexp_topN(set1_filter, set2_filter, count)
        return make_response(diffexp, HTTPStatus.OK, {"Content-Type": "application/json"})
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)
    except JSONEncodingValueError:
        # JSON encoding failure, usually due to bad data. Just let it ripple up
        # to default exception handler.
        current_app.logger.warning(JSON_NaN_to_num_warning_msg)
        raise


def layout_obs_get(request, data_adaptor):
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    try:
        if preferred_mimetype == "application/octet-stream":
            return make_response(
                data_adaptor.layout_to_fbs_matrix(), HTTPStatus.OK, {"Content-Type": "application/octet-stream"}
            )
        else:
            return abort(HTTPStatus.NOT_ACCEPTABLE)
    except PrepareError:
        return abort_and_log(
            HTTPStatus.NOT_IMPLEMENTED,
            f"No embedding available {request.path}",
            loglevel=logging.ERROR,
            include_exc_info=True,
        )


def layout_obs_put(request, data_adaptor):
    if not data_adaptor.config.enable_reembedding:
        return abort(HTTPStatus.NOT_IMPLEMENTED)

    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    args = request.get_json()
    filter = args["filter"] if args else None
    if not filter:
        return abort_and_log(HTTPStatus.BAD_REQUEST, "obs filter is required")
    method = args["method"] if args else "umap"

    try:
        schema, fbs = data_adaptor.compute_embedding(method, filter)
        return make_response(
            fbs,
            HTTPStatus.OK,
            {
                "Content-Type": "application/octet-stream",
                "CxG-Schema": json.dumps(schema),
                "Access-Control-Expose-Headers": "CxG-Schema",
            },
        )
    except NotImplementedError as e:
        return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, str(e), include_exc_info=True)
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)
