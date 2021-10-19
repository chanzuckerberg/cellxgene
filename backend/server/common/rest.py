import copy
import logging
import sys
from http import HTTPStatus
import zlib
import json
import pandas as pd
import numpy as np
from flask import make_response, jsonify, current_app, abort
from werkzeug.urls import url_unquote

from backend.server.common.config.client_config import get_client_config, get_client_userinfo
from backend.common.constants import Axis, DiffExpMode, JSON_NaN_to_num_warning_msg
from backend.common.errors import (
    FilterError,
    JSONEncodingValueError,
    PrepareError,
    DisabledFeatureError,
    ExceedsLimitError,
    DatasetAccessError,
    ColorFormatException,
    AnnotationsError,
    ObsoleteRequest,
    UnsupportedSummaryMethod,
)
from backend.common.genesets import summarizeQueryHash
from backend.common.fbs.matrix import decode_matrix_fbs


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


def _query_parameter_to_filter(args):
    """
    Convert an annotation value filter, if present in the query args,
    into the standard dict filter format used by internal code.

    Query param filters look like:  <axis>:name=value, where value
    may be one of:
        - a range, min,max, where either may be an open range by using an asterisk, eg, 10,*
        - a value
    Eg,
        ...?tissue=lung&obs:tissue=heart&obs:num_reads=1000,*
    """
    filters = {
        "obs": {},
        "var": {},
    }

    # args has already been url-unquoted once.  We assume double escaping
    # on name and value.
    try:
        for key, value in args.items(multi=True):
            axis, name = key.split(":")
            if axis not in ("obs", "var"):
                raise FilterError("unknown filter axis")
            name = url_unquote(name)
            current = filters[axis].setdefault(name, {"name": name})

            val_split = value.split(",")
            if len(val_split) == 1:
                if "min" in current or "max" in current:
                    raise FilterError("do not mix range and value filters")
                value = url_unquote(value)
                values = current.setdefault("values", [])
                values.append(value)

            elif len(val_split) == 2:
                if len(current) > 1:
                    raise FilterError("duplicate range specification")
                min = url_unquote(val_split[0])
                max = url_unquote(val_split[1])
                if min != "*":
                    current["min"] = float(min)
                if max != "*":
                    current["max"] = float(max)
                if len(current) < 2:
                    raise FilterError("must specify at least min or max in range filter")

            else:
                raise FilterError("badly formated filter value")

    except ValueError as e:
        raise FilterError(str(e))

    result = {}
    for axis in ("obs", "var"):
        axis_filter = filters[axis]
        if len(axis_filter) > 0:
            result[axis] = {"annotation_value": [val for val in axis_filter.values()]}

    return result


def schema_get_helper(data_adaptor):
    """helper function to gather the schema from the data source and annotations"""
    schema = data_adaptor.get_schema()
    schema = copy.deepcopy(schema)

    # add label obs annotations as needed
    annotations = data_adaptor.dataset_config.user_annotations
    if annotations.user_annotations_enabled():
        label_schema = annotations.get_schema(data_adaptor)
        curr_schema = schema["annotations"]["obs"]["columns"]
        keys = []
        vals = []
        for i in (label_schema + curr_schema):
            if i['name'] not in keys:
                keys.append(i['name'])
                vals.append(i)
        schema["annotations"]["obs"]["columns"] = vals

    return schema


def schema_get(data_adaptor):
    schema = schema_get_helper(data_adaptor)
    return make_response(jsonify({"schema": schema}), HTTPStatus.OK)


def config_get(app_config, data_adaptor):
    config = get_client_config(app_config, data_adaptor)
    return make_response(jsonify(config), HTTPStatus.OK)


def userinfo_get(app_config, data_adaptor):
    config = get_client_userinfo(app_config, data_adaptor)
    return make_response(jsonify(config), HTTPStatus.OK)


def annotations_obs_get(request, data_adaptor):
    fields = request.args.getlist("annotation-name", None)
    num_columns_requested = len(data_adaptor.get_obs_keys()) if len(fields) == 0 else len(fields)
    if data_adaptor.server_config.exceeds_limit("column_request_max", num_columns_requested):
        return abort(HTTPStatus.BAD_REQUEST)
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    try:
        labels = None
        labels = data_adaptor.data.obs[fields]
        fbs = data_adaptor.annotation_to_fbs_matrix(Axis.OBS, fields, labels)
        return make_response(fbs, HTTPStatus.OK, {"Content-Type": "application/octet-stream"})
    except KeyError as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def annotations_put_fbs_helper(data_adaptor, fbs):
    """helper function to write annotations from fbs"""
    annotations = data_adaptor.dataset_config.user_annotations
    if not annotations.user_annotations_enabled():
        raise DisabledFeatureError("Writable annotations are not enabled")

    new_label_df = decode_matrix_fbs(fbs)
    if not new_label_df.empty:
        new_label_df = data_adaptor.check_new_labels(new_label_df)

    for k in new_label_df.columns:
        data_adaptor.data.obs[k] = pd.Categorical(np.array(list(new_label_df[k])).astype('str'))
    
    data_adaptor._save_orig_data()


def inflate(data):
    return zlib.decompress(data)

def annotations_obs_put(request, data_adaptor):
    annotations = data_adaptor.dataset_config.user_annotations
    if not annotations.user_annotations_enabled():
        return abort(HTTPStatus.NOT_IMPLEMENTED)

    anno_collection = request.args.get("annotation-collection-name", default=None)
    fbs = inflate(request.get_data())

    if anno_collection is not None:
        if not annotations.is_safe_collection_name(anno_collection):
            return abort(HTTPStatus.BAD_REQUEST, "Bad annotation collection name")
        annotations.set_collection(anno_collection)

    try:
        annotations_put_fbs_helper(data_adaptor, fbs)
        res = json.dumps({"status": "OK"})
        data_adaptor._create_schema()
        return make_response(res, HTTPStatus.OK, {"Content-Type": "application/json"})
    except (ValueError, DisabledFeatureError, KeyError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def annotations_var_get(request, data_adaptor):
    fields = request.args.getlist("annotation-name", None)
    num_columns_requested = len(data_adaptor.get_var_keys()) if len(fields) == 0 else len(fields)
    if data_adaptor.server_config.exceeds_limit("column_request_max", num_columns_requested):
        return abort(HTTPStatus.BAD_REQUEST)
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    try:
        labels = None
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
    except (FilterError, ValueError, ExceedsLimitError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def data_var_get(request, data_adaptor):
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    try:
        filter = _query_parameter_to_filter(request.args)
        return make_response(
            data_adaptor.data_frame_to_fbs_matrix(filter, axis=Axis.VAR),
            HTTPStatus.OK,
            {"Content-Type": "application/octet-stream"},
        )
    except (FilterError, ValueError, ExceedsLimitError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def colors_get(data_adaptor):
    if not data_adaptor.dataset_config.presentation__custom_colors:
        return make_response(jsonify({}), HTTPStatus.OK)
    try:
        return make_response(jsonify(data_adaptor.get_colors()), HTTPStatus.OK)
    except ColorFormatException as e:
        return abort_and_log(HTTPStatus.NOT_FOUND, str(e), include_exc_info=True)


def diffexp_obs_post(request, data_adaptor):
    if not data_adaptor.dataset_config.diffexp__enable:
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
    except (ValueError, DisabledFeatureError, FilterError, ExceedsLimitError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)
    except JSONEncodingValueError:
        # JSON encoding failure, usually due to bad data. Just let it ripple up
        # to default exception handler.
        current_app.logger.warning(JSON_NaN_to_num_warning_msg)
        raise


def layout_obs_get(request, data_adaptor):
    fields = request.args.getlist("layout-name", None)
    num_columns_requested = len(data_adaptor.get_embedding_names()) if len(fields) == 0 else len(fields)
    if data_adaptor.server_config.exceeds_limit("column_request_max", num_columns_requested):
        return abort(HTTPStatus.BAD_REQUEST)

    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    try:
        return make_response(
            data_adaptor.layout_to_fbs_matrix(fields), HTTPStatus.OK, {"Content-Type": "application/octet-stream"}
        )
    except (KeyError, DatasetAccessError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)
    except PrepareError:
        return abort_and_log(
            HTTPStatus.NOT_IMPLEMENTED,
            f"No embedding available {request.path}",
            loglevel=logging.ERROR,
            include_exc_info=True,
        )

def sankey_data_put(request, data_adaptor):
    args = request.get_json()
    labels = args.get("labels", None)
    labels = [x['__columns'][0] for x in labels]
    name = args.get("name", None)

    if not labels:
        return make_response(jsonify({"edges":[],"weights":[]}),HTTPStatus.OK, {"Content-Type": "application/json"})
    

    try:
        edges,weights = data_adaptor.compute_sankey_df(labels,name)
        edges = [list(x) for x in edges]
        weights = list(weights)
        return make_response(jsonify({"edges": edges, "weights": weights}), HTTPStatus.OK, {"Content-Type": "application/json"})
    except NotImplementedError as e:
        return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, str(e))
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)

def output_data_put(request, data_adaptor):
    args = request.get_json()
    saveName = args.get("saveName","")
    labels = args.get("labels",None)
    labelNames = args.get("labelNames",None)
    currentLayout = args.get("currentLayout",None)
    filter = args["filter"] if args else None
    
    try:
        shape = data_adaptor.get_shape()
        obs_mask = data_adaptor._axis_filter_to_mask(Axis.OBS, filter["obs"], shape[0])
    except (KeyError, IndexError):
        raise FilterError("Error parsing filter")

    fail=False
    if saveName != "":
        X = data_adaptor.data.obsm["X_"+currentLayout]
        f = np.isnan(X).sum(1)==0    
        filt = np.logical_and(f,obs_mask)

        temp = {}
        for key in data_adaptor.data.uns.keys():
            if key[:2] == "N_" and "_mask" not in key:
                temp[key] = data_adaptor.data.uns[key][filt][:,filt]
            elif key[:2] == "N_":
                temp[key] = data_adaptor.data.uns[key][filt]

        adata = data_adaptor.data[filt].copy()

        for key in temp.keys():
            adata.uns[key] = temp[key]
            
        name = currentLayout.split(';')[-1]
        
        if labels and labelNames:
            labels = [x['__columns'][0] for x in labels]
            for n,l in zip(labelNames,labels):
                adata.obs[n] = pd.Categorical(l)        
        
        keys = list(adata.obsm.keys())
        for k in keys:
            if k[:2] == "X_":
                k2 = k[2:]
            else:
                k2 = k
            if name not in k2.split(';;'):
                del adata.obsm[k]
        
        keys = list(adata.uns.keys())
        for k in keys:
            if k[:2] == "N_":
                k2 = k[2:]
            else:
                k2 = k
            k2 = k2.split('_mask')[0]
            if name not in k2.split(';;') and k[:2] == "N_":          
                del adata.uns[k]  

        try:
            adata.obs_names = pd.Index(adata.obs["name_0"].astype('str'))
            del adata.obs["name_0"]
        except:
            pass
        try:
            adata.var_names = pd.Index(adata.var["name_0"].astype('str'))
            del adata.var["name_0"]
        except:
            pass        

        if "X" in adata.layers.keys():
            adata.X = adata.layers["X"]
            del adata.layers["X"]
        
        try:
            del adata.obsm["X_root"]
        except:
            pass

        try:
            adata.write_h5ad(saveName.split('.h5ad')[0]+'.h5ad')
        except:
            fail=True

    try:
        return make_response(jsonify({"fail": fail}), HTTPStatus.OK, {"Content-Type": "application/json"})
    except NotImplementedError as e:
        return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, str(e))
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)  

def reload_put(request, data_adaptor):
    args = request.get_json()
    currentLayout = args.get("currentLayout",None)
    filter = args["filter"] if args else None
    
    try:
        shape = data_adaptor.get_shape()
        obs_mask = data_adaptor._axis_filter_to_mask(Axis.OBS, filter["obs"], shape[0])
    except (KeyError, IndexError):
        raise FilterError("Error parsing filter")

    fail=False
    
    X = data_adaptor.data.obsm["X_"+currentLayout]
    f = np.isnan(X).sum(1)==0    
    filt = np.logical_and(f,obs_mask)

    temp = {}
    for key in data_adaptor.data.uns.keys():
        if key[:2] == "N_" and "_mask" not in key:
            temp[key] = data_adaptor.data.uns[key][filt][:,filt]
        elif key[:2] == "N_":
            temp[key] = data_adaptor.data.uns[key][filt]

    adata = data_adaptor.data[filt].copy()

    for key in temp.keys():
        adata.uns[key] = temp[key]

    data_adaptor.data = adata
    
    try:
        data_adaptor.data.obs_names = pd.Index(data_adaptor.data.obs["name_0"].astype('str'))
        del data_adaptor.data.obs["name_0"]
    except:
        pass
    try:
        data_adaptor.data.var_names = pd.Index(data_adaptor.data.var["name_0"].astype('str'))
        del data_adaptor.data.var["name_0"]
    except:
        pass        

    data_adaptor._validate_and_initialize()
    
    try:
        return make_response(jsonify({"fail": fail}), HTTPStatus.OK, {"Content-Type": "application/json"})
    except NotImplementedError as e:
        return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, str(e))
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)  

def reload_full_put(request, data_adaptor):
    data_adaptor.data = data_adaptor.data_orig
    try:
        data_adaptor.data.obs_names = pd.Index(data_adaptor.data.obs["name_0"].astype('str'))
        del data_adaptor.data.obs["name_0"]
    except:
        pass
    try:
        data_adaptor.data.var_names = pd.Index(data_adaptor.data.var["name_0"].astype('str'))
        del data_adaptor.data.var["name_0"]
    except:
        pass        

    data_adaptor._validate_and_initialize()
    data_adaptor.data_orig = data_adaptor.data
    
    try:
        return make_response(jsonify({"fail": False}), HTTPStatus.OK, {"Content-Type": "application/json"})
    except NotImplementedError as e:
        return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, str(e))
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)  

def rename_obs_put(request, data_adaptor):
    args = request.get_json()
    oldName = args.get("oldCategoryName","")
    newName = args.get("newCategoryName","")
    

    fail=False
    if oldName != "" and newName != "":
        try:
            data_adaptor.data.obs[newName] = data_adaptor.data.obs[oldName]
            del data_adaptor.data.obs[oldName]         
            data_adaptor._create_schema()
        except:
            fail=True
    
    data_adaptor._save_orig_data()
    try:
        return make_response(jsonify({"fail": fail}), HTTPStatus.OK, {"Content-Type": "application/json"})
    except NotImplementedError as e:
        return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, str(e))
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)  

def delete_obs_put(request, data_adaptor):
    args = request.get_json()
    name = args.get("category","")

    fail=False
    if name != "":
        try:
            del data_adaptor.data.obs[name]
            data_adaptor._create_schema()
        except:
            fail=True
    
    data_adaptor._save_orig_data()
    try:
        return make_response(jsonify({"fail": fail}), HTTPStatus.OK, {"Content-Type": "application/json"})
    except NotImplementedError as e:
        return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, str(e))
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)  

def delete_obsm_put(request, data_adaptor):
    args = request.get_json()
    embNames = args.get("embNames",None)
    fail=False
    if embNames is not None:
        for embName in embNames:
            try:
                del data_adaptor.data.obsm[f'X_{embName}']
            except:
                pass
            try:
                del data_adaptor.data.uns[f"N_{embName}"]
                del data_adaptor.data.uns[f"N_{embName}_mask"]
            except:
                pass

            data_adaptor._refresh_layout_schema()
    
    data_adaptor._save_orig_data()            
    try:
        return make_response(jsonify({"fail": fail}), HTTPStatus.OK, {"Content-Type": "application/json"})
    except NotImplementedError as e:
        return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, str(e))
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)    

def rename_obsm_put(request, data_adaptor):
    args = request.get_json()
    embNames = args.get("embNames",None)
    oldName = args.get("oldName",None)
    newName = args.get("newName",None)

    if embNames is not None and oldName is not None and newName is not None:
        for embName in embNames:
            try:
                X1 = data_adaptor.data.obsm[f'X_{embName}']
                del data_adaptor.data.obsm[f'X_{embName}']
                data_adaptor.data.obsm[f'X_{embName.replace(oldName,newName)}'] = X1
            except:
                pass
            try:
                X1 = data_adaptor.data.uns[f"N_{embName}"]
                X2 = data_adaptor.data.uns[f"N_{embName}_mask"]
                del data_adaptor.data.uns[f"N_{embName}"]
                del data_adaptor.data.uns[f"N_{embName}_mask"]
                data_adaptor.data.uns[f"N_{embName.replace(oldName,newName)}"] = X1
                data_adaptor.data.uns[f"N_{embName.replace(oldName,newName)}_mask"] = X2
            except:
                pass

            data_adaptor._refresh_layout_schema()
    
    data_adaptor._save_orig_data()            
    try:
        return make_response(jsonify({"schema": data_adaptor.get_schema()}), HTTPStatus.OK, {"Content-Type": "application/json"})
    except NotImplementedError as e:
        return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, str(e))
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)    


def change_layer_put(request, data_adaptor):
    args = request.get_json()
    dataLayer = args.get("dataLayer",None)
    fail=False
    if dataLayer is not None:            
        try:
            if "X" not in data_adaptor.data.layers.keys():
                data_adaptor.data.layers['X'] = data_adaptor.data.X
            if dataLayer == ".raw":
                data_adaptor.data.X = data_adaptor.data.raw.X
            else:
                data_adaptor.data.X = data_adaptor.data.layers[dataLayer]
        except:
            fail=True
    try:
        return make_response(jsonify({"fail": fail}), HTTPStatus.OK, {"Content-Type": "application/json"})
    except NotImplementedError as e:
        return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, str(e))
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)    

def leiden_put(request, data_adaptor):
    args = request.get_json()
    name = args.get("name", None)
    cName = args.get("cName", None)
    resolution = args.get('resolution')

    try:
        cl = list(data_adaptor.compute_leiden(name,cName,resolution))
        return make_response(jsonify({"clusters": cl}), HTTPStatus.OK, {"Content-Type": "application/json"})
    except NotImplementedError as e:
        return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, str(e))
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def layout_obs_put(request, data_adaptor):
    args = request.get_json()
    filter = args["filter"] if args else None
    if not filter:
        return abort_and_log(HTTPStatus.BAD_REQUEST, "obs filter is required")
    method = args["method"] if args else "umap"
    reembedParams = args["params"] if args else {}
    parentName = args["parentName"] if args else ""
    embName = args["embName"] if args else None

    try:
        schema = data_adaptor.compute_embedding(method, filter, reembedParams, parentName, embName)
        return make_response(jsonify(schema), HTTPStatus.OK, {"Content-Type": "application/json"})
    except NotImplementedError as e:
        return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, str(e))
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def preprocess_put(request, data_adaptor):
    args = request.get_json()
    reembedParams = args["params"] if args else {}

    try:
        schema = data_adaptor.compute_preprocess(reembedParams)
        return make_response(jsonify(schema), HTTPStatus.OK, {"Content-Type": "application/json"})
    except NotImplementedError as e:
        return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, str(e))
    except (ValueError, DisabledFeatureError, FilterError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def reembed_parameters_get(request, data_adaptor):
    preferred_mimetype = request.accept_mimetypes.best_match(["application/json", "text/csv"])
    if preferred_mimetype not in ("application/json", "text/csv"):
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    try:
        annotations = data_adaptor.dataset_config.user_annotations
        (reembedParams) = annotations.read_reembed_parameters(data_adaptor)

        return make_response(
            jsonify({"reembedParams": reembedParams}), HTTPStatus.OK
        )
    except (ValueError, KeyError, AnnotationsError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e))

def genesets_get(request, data_adaptor):
    preferred_mimetype = request.accept_mimetypes.best_match(["application/json", "text/csv"])
    if preferred_mimetype not in ("application/json", "text/csv"):
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    try:
        annotations = data_adaptor.dataset_config.user_annotations
        (genesets, tid) = annotations.read_gene_sets(data_adaptor)

        if preferred_mimetype == "text/csv":
            return make_response(
                annotations.gene_sets_to_csv(genesets),
                HTTPStatus.OK,
                {
                    "Content-Type": "text/csv",
                    "Content-Disposition": "attachment; filename=genesets.csv",
                },
            )
        else:
            return make_response(
                jsonify({"genesets": annotations.gene_sets_to_response(genesets), "tid": tid}), HTTPStatus.OK
            )
    except (ValueError, KeyError, AnnotationsError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e))

def reembed_parameters_put(request, data_adaptor):
    annotations = data_adaptor.dataset_config.user_annotations

    anno_collection = request.args.get("annotation-collection-name", default=None)
    if anno_collection is not None:
        if not annotations.is_safe_collection_name(anno_collection):
            return abort(HTTPStatus.BAD_REQUEST, "Bad annotation collection name")
        annotations.set_collection(anno_collection)

    args = request.get_json()
    try:
        reembedParams = args.get("reembedParams", None)
        if reembedParams is None:
            abort(HTTPStatus.BAD_REQUEST)

        annotations.write_reembed_parameters(reembedParams, data_adaptor)
        return make_response(jsonify({"status": "OK"}), HTTPStatus.OK)
    except (ValueError, DisabledFeatureError, KeyError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)
    except (ObsoleteRequest, TypeError) as e:
        return abort(HTTPStatus.NOT_FOUND, description=str(e))


def genesets_put(request, data_adaptor):
    annotations = data_adaptor.dataset_config.user_annotations
    if not annotations.gene_sets_save_enabled():
        return abort(HTTPStatus.NOT_IMPLEMENTED)

    anno_collection = request.args.get("annotation-collection-name", default=None)
    if anno_collection is not None:
        if not annotations.is_safe_collection_name(anno_collection):
            return abort(HTTPStatus.BAD_REQUEST, "Bad annotation collection name")
        annotations.set_collection(anno_collection)

    args = request.get_json()
    try:
        genesets = args.get("genesets", None)
        tid = args.get("tid", None)
        if genesets is None:
            abort(HTTPStatus.BAD_REQUEST)

        annotations.write_gene_sets(genesets, tid, data_adaptor)
        return make_response(jsonify({"status": "OK"}), HTTPStatus.OK)
    except (ValueError, DisabledFeatureError, KeyError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)
    except (ObsoleteRequest, TypeError) as e:
        return abort(HTTPStatus.NOT_FOUND, description=str(e))


def summarize_var_helper(request, data_adaptor, key, raw_query):
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    summary_method = request.values.get("method", default="mean")
    query_hash = summarizeQueryHash(raw_query)
    if key and query_hash != key:
        return abort(HTTPStatus.BAD_REQUEST, description="query key did not match")

    args_filter_only = request.values.copy()
    args_filter_only.poplist("method")
    args_filter_only.poplist("key")

    try:
        filter = _query_parameter_to_filter(args_filter_only)
        return make_response(
            data_adaptor.summarize_var(summary_method, filter, query_hash),
            HTTPStatus.OK,
            {"Content-Type": "application/octet-stream"},
        )
    except (ValueError) as e:
        return abort(HTTPStatus.NOT_FOUND, description=str(e))
    except (UnsupportedSummaryMethod, FilterError) as e:
        return abort(HTTPStatus.BAD_REQUEST, description=str(e))


def summarize_var_get(request, data_adaptor):
    return summarize_var_helper(request, data_adaptor, None, request.query_string)


def summarize_var_post(request, data_adaptor):
    if not request.content_type or "application/x-www-form-urlencoded" not in request.content_type:
        return abort(HTTPStatus.UNSUPPORTED_MEDIA_TYPE)
    if request.content_length > 1_000_000:  # just a sanity check to avoid memory exhaustion
        return abort(HTTPStatus.BAD_REQUEST)

    key = request.args.get("key", default=None)
    return summarize_var_helper(request, data_adaptor, key, request.get_data())
