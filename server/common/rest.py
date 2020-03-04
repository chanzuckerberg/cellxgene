from http import HTTPStatus
import warnings
import copy
from flask import make_response, jsonify
from server.common.constants import Axis, DiffExpMode, JSON_NaN_to_num_warning_msg
from server.common.errors import (
    FilterError,
    JSONEncodingValueError,
    PrepareError,
    DisabledFeatureError,
)

import json
from server.data_common.fbs.matrix import decode_matrix_fbs


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
        return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
    try:
        labels = None
        if annotations:
            labels = annotations.read_labels(data_adaptor)
        fbs = data_adaptor.annotation_to_fbs_matrix(Axis.OBS, fields, labels)
        return make_response(fbs, HTTPStatus.OK, {"Content-Type": "application/octet-stream"})
    except KeyError:
        return make_response(f"Error bad key in {fields}", HTTPStatus.BAD_REQUEST)
    except ValueError as e:
        return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)
    except Exception as e:
        return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


def annotations_put_fbs_helper(data_adaptor, annotations, fbs):
    """helper function to write annotations from fbs"""
    if annotations is None:
        raise DisabledFeatureError("Writable annotations are not enabled")

    new_label_df = decode_matrix_fbs(fbs)
    if not new_label_df.empty:
        data_adaptor.check_new_labels(new_label_df)
    annotations.write_labels(new_label_df, data_adaptor)


def annotations_obs_put(request, data_adaptor, annotations):
    anno_collection = request.args.get("annotation-collection-name", default=None)
    fbs = request.get_data()
    if annotations is None:
        return make_response("Error, annotations are not configured", HTTPStatus.BAD_REQUEST)

    if anno_collection is not None:
        if not annotations.is_safe_collection_name(anno_collection):
            return make_response(f"Error, bad annotation collection name", HTTPStatus.BAD_REQUEST)
        annotations.set_collection(anno_collection)

    try:
        annotations_put_fbs_helper(data_adaptor, annotations, fbs)
        res = json.dumps({"status": "OK"})
        return make_response(res, HTTPStatus.OK, {"Content-Type": "application/json"})
    except (ValueError, DisabledFeatureError, KeyError) as e:
        return make_response(str(e), HTTPStatus.BAD_REQUEST)
    except Exception as e:
        return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


def annotations_var_get(request, data_adaptor, annotations):
    fields = request.args.getlist("annotation-name", None)
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
    try:
        labels = None
        if annotations is not None:
            labels = annotations.read_labels(data_adaptor)
        return make_response(
            data_adaptor.annotation_to_fbs_matrix(Axis.VAR, fields, labels),
            HTTPStatus.OK,
            {"Content-Type": "application/octet-stream"},
        )
    except KeyError:
        return make_response(f"Error bad key in {fields}", HTTPStatus.BAD_REQUEST)
    except ValueError as e:
        return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)
    except Exception as e:
        return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


def data_var_put(request, data_adaptor):
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)

    filter_json = request.get_json()
    filter = filter_json["filter"] if filter_json else None
    try:
        return make_response(
            data_adaptor.data_frame_to_fbs_matrix(filter, axis=Axis.VAR),
            HTTPStatus.OK,
            {"Content-Type": "application/octet-stream"},
        )
    except FilterError as e:
        return make_response(str(e), HTTPStatus.BAD_REQUEST)
    except ValueError as e:
        return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


def diffexp_obs_post(request, data_adaptor):
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
        diffexp = data_adaptor.diffexp_topN(set1_filter, set2_filter, count)
        return make_response(diffexp, HTTPStatus.OK, {"Content-Type": "application/json"})
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return make_response(str(e), HTTPStatus.BAD_REQUEST)
    except JSONEncodingValueError as e:
        # JSON encoding failure, usually due to bad data
        warnings.warn(JSON_NaN_to_num_warning_msg)
        return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)
    except ValueError as e:
        return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


def layout_obs_get(request, data_adaptor):
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    try:
        if preferred_mimetype == "application/octet-stream":
            return make_response(
                data_adaptor.layout_to_fbs_matrix(), HTTPStatus.OK, {"Content-Type": "application/octet-stream"}
            )
        else:
            return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
    except PrepareError as e:
        return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)
    except ValueError as e:
        return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


def layout_obs_put(request, data_adaptor):
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
    if not data_adaptor.config.enable_reembedding:
        return make_response(f"Computed embedding not supported.", HTTPStatus.BAD_REQUEST)

    args = request.get_json()
    filter = args["filter"] if args else None
    if not filter:
        return make_response("Error: obs filter is required", HTTPStatus.BAD_REQUEST)
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
        return make_response(str(e), HTTPStatus.NOT_IMPLEMENTED)
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return make_response(str(e), HTTPStatus.BAD_REQUEST)
