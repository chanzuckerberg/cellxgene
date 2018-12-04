from http import HTTPStatus
import pkg_resources
import warnings

from flask import (
    Blueprint, current_app, jsonify, make_response, request
)
from flask_restful_swagger_2 import Api, swagger, Resource
from werkzeug.datastructures import ImmutableMultiDict

from server.app.util.constants import Axis, DiffExpMode, JSON_NaN_to_num_warning_msg
from server.app.util.filter import parse_filter, QueryStringError
from server.app.util.models import FilterModel
from server.app.util.utils import get_mime_type
from server.app.util.errors import MimeTypeError, FilterError, InteractiveError, PrepareError

"""
Sort order for routes
1. Initialize
2. Data & Metadata
3. Computation
"""


class SchemaAPI(Resource):
    @swagger.doc({
        "summary": "get schema for dataframe and annotations",
        "tags": ["initialize"],
        "parameters": [],
        "responses": {
            "200": {
                "description": "schema",
                "examples": {
                    "application/json": {
                        "schema": {
                            "dataframe": {
                                "nObs": 383,
                                "nVar": 19944,
                                "type": "float32"
                            },
                            "annotations": {
                                "obs": [
                                    {"name": "name", "type": "string"},
                                    {"name": "tissue_type", "type": "string"},
                                    {"name": "num_reads", "type": "int32"},
                                    {"name": "sample_name", "type": "string"},
                                    {
                                        "name": "clusters",
                                        "type": "categorical",
                                        "categories": [99, 1, "unknown cluster"]
                                    },
                                    {"name": "QScore", "type": "float32"}
                                ],
                                "var": [
                                    {"name": "name", "type": "string"},
                                    {"name": "gene", "type": "string"}
                                ]
                            }
                        }
                    }
                }
            }
        }

    })
    def get(self):
        return make_response(jsonify({"schema": current_app.data.schema}), HTTPStatus.OK)


class ConfigAPI(Resource):
    @swagger.doc({
        "summary": "Configuration information to assist in front-end adaptation"
                   " to underlying engine, available functionality, interactive time limits, etc",
        "tags": ["initialize"],
        "parameters": [],
        "responses": {
            "200": {
                "description": "schema",
                "examples": {
                    "application/json": {
                        "config": {
                            "features": [
                                {"method": "POST", "path": "/cluster/", "available": False},
                                {
                                    "method": "POST",
                                    "path": "/layout/obs",
                                    "available": True,
                                    "interactiveLimit": 10000
                                },
                                {"method": "POST", "path": "/layout/var", "available": False}

                            ],
                            "displayNames": {
                                "engine": "ScanPy version 1.33",
                                "dataset": "/home/joe/mouse/blorth.csv"
                            },
                        }
                    }
                }
            }
        }
    })
    def get(self):
        config = {
            "config": {
                "features": [
                    {"method": "POST", "path": "/cluster/", **current_app.data.features["cluster"]},
                    {"method": "POST", "path": "/layout/obs", **current_app.data.features["layout"]["obs"]},
                    {"method": "POST", "path": "/layout/var", **current_app.data.features["layout"]["var"]},
                    {"method": "POST", "path": "/diffexp/", **current_app.data.features["diffexp"]},
                ],
                "displayNames": {
                    "engine": f"cellxgene Scanpy engine version {pkg_resources.get_distribution('cellxgene').version}",
                    "dataset": current_app.config["DATASET_TITLE"]
                },
                "parameters": {
                    "max_category_items": current_app.data.max_category_items
                }
            }
        }
        return make_response(jsonify(config), HTTPStatus.OK)


class AnnotationsObsAPI(Resource):
    @swagger.doc({
        "summary": "Fetch annotations (metadata) for all observations.",
        "tags": ["annotations"],
        "parameters": [{
            "in": "query",
            "name": "annotation-name",
            "type": "string",
            "description": "list of 1 or more annotation names"
        }],
        "responses": {
            "200": {
                "description": "annotations",
                "examples": {
                    "application/json": {
                        "names": [
                            "tissue_type", "sex", "num_reads", "clusters"
                        ],
                        "data": [
                            [0, "lung", "F", 39844, 99],
                            [1, "heart", "M", 83, 1],
                            [49, "spleen", None, 2, "unknown cluster"],

                        ]
                    }

                }
            },
            "400": {
                "description": "one or more of the annotation-name identifiers were not associated with an "
                               "annotation name"
            }
        }
    })
    def get(self):
        fields = request.args.getlist("annotation-name", None)
        try:
            annotation_response = current_app.data.annotation({}, "obs", fields)
        except KeyError:
            return make_response(f"Error bad key in {fields}", HTTPStatus.BAD_REQUEST)
        try:
            return make_response(jsonify(annotation_response), HTTPStatus.OK)
        except ValueError as e:
            # JSON encoding failure, usually due to bad data
            warnings.warn(JSON_NaN_to_num_warning_msg)
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)

    @swagger.doc({
        "summary": "Fetch annotations (metadata) for filtered subset of observations.",
        "tags": ["annotations"],
        "parameters": [
            {
                "in": "query",
                "name": "annotation-name",
                "type": "string",
                "description": "list of 1 or more annotation names"
            },
            {
                "name": "filter",
                "description": "Complex Filter",
                "in": "body",
                "schema": FilterModel
            }
        ],
        "responses": {
            "200": {
                "description": "annotations",
                "examples": {
                    "application/json": {
                        "names": [
                            "tissue_type", "sex", "num_reads", "clusters"
                        ],
                        "data": [
                            [0, "lung", "F", 39844, 99],
                            [1, "heart", "M", 83, 1],
                            [49, "spleen", None, 2, "unknown cluster"],

                        ]
                    }

                }
            },
            "400": {
                "description": "malformed filter or one or more of the annotation-name identifiers were"
                               "not associated with an annotation name"
            }
        }
    })
    def put(self):
        fields = request.args.getlist("annotation-name", None)
        try:
            annotation_response = current_app.data.annotation(request.get_json()["filter"], "obs", fields)
        except KeyError:
            return make_response(f"Error bad key in {fields}", HTTPStatus.BAD_REQUEST)
        except FilterError as e:
            return make_response(e.message, HTTPStatus.BAD_REQUEST)
        try:
            return make_response(jsonify(annotation_response), HTTPStatus.OK)
        except ValueError as e:
            # JSON encoding failure, usually due to bad data
            warnings.warn(JSON_NaN_to_num_warning_msg)
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


class AnnotationsVarAPI(Resource):
    @swagger.doc({
        "summary": "Fetch annotations (metadata) for all variables.",
        "tags": ["annotations"],
        "parameters": [{
            "in": "query",
            "name": "annotation-name",
            "type": "string",
            "description": "list of 1 or more annotation names"
        }],
        "responses": {
            "200": {
                "description": "annotations",
                "examples": {
                    "application/json": {
                        "names": [
                            "name", "category"
                        ],
                        "data": [
                            [0, "ATAD3C", 1],
                            [1, "RER1", None],
                            [49, "S100B", 6]
                        ]
                    }

                }
            },
            "400": {
                "description": "one or more of the annotation-name identifiers were not associated with an"
                               " annotation name"
            }
        }
    })
    def get(self):
        fields = request.args.getlist("annotation-name", None)
        try:
            annotation_response = current_app.data.annotation({}, "var", fields)
        except KeyError:
            return make_response(f"Error bad key in {fields}", HTTPStatus.BAD_REQUEST)
        try:
            return make_response(jsonify(annotation_response), HTTPStatus.OK)
        except ValueError as e:
            # JSON encoding failure, usually due to bad data
            warnings.warn(JSON_NaN_to_num_warning_msg)
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)

    @swagger.doc({
        "summary": "Fetch annotations (metadata) for filtered subset of variables.",
        "tags": ["annotations"],
        "parameters": [
            {
                "in": "query",
                "name": "annotation-name",
                "type": "string",
                "description": "list of 1 or more annotation names"
            },
            {
                "name": "filter",
                "description": "Complex Filter",
                "in": "body",
                "schema": FilterModel
            }
        ],
        "responses": {
            "200": {
                "description": "annotations",
                "examples": {
                    "application/json": {
                        "names": [
                            "name", "category"
                        ],
                        "data": [
                            [0, "ATAD3C", 1],
                            [1, "RER1", None],
                            [49, "S100B", 6]
                        ]
                    }
                }
            },
            "400": {
                "description": "malformed filter or one or more of the annotation-name identifiers were"
                               "not associated with an annotation name"
            }
        }
    })
    def put(self):
        fields = request.args.getlist("annotation-name", None)
        try:
            annotation_response = current_app.data.annotation(request.get_json()["filter"], "var", fields)
        except KeyError:
            return make_response(f"Error bad key in {fields}", HTTPStatus.BAD_REQUEST)
        except FilterError:
            return make_response("Malformed filter", HTTPStatus.BAD_REQUEST)
        try:
            return make_response(jsonify(annotation_response), HTTPStatus.OK)
        except ValueError as e:
            # JSON encoding failure, usually due to bad data
            warnings.warn(JSON_NaN_to_num_warning_msg)
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


class DataObsAPI(Resource):
    @swagger.doc({
        "summary": "Get data (expression values) from the dataframe.",
        "tags": ["data"],
        "parameters": [
            {
                "in": "query",
                "name": "filter",
                "type": "string",
                "description": "axis:key:value"
            },
            {
                "in": "query",
                "name": "accept-type",
                "type": "string",
                "description": "MIME type"
            },
        ],
        "responses": {
            "200": {
                "description": "expression",
                "examples": {
                    "application/json": {
                        "var": [0, 20000],
                        "obs": [
                            [1, 39483, 3902, 203, 0, 0, 28]
                        ]
                    }
                }
            },
            "400": {
                "description": "Malformed filter"
            },
            "406": {
                "description": "Unacceptable MIME type"
            },
        }
    })
    def get(self):
        accept_type = request.args.get("accept-type", None)
        # request.args is immutable
        args = request.args.copy()
        args.pop("accept-type", None)
        try:
            filter_ = parse_filter(ImmutableMultiDict(args), current_app.data.schema['annotations'])
        except QueryStringError as e:
            return make_response(e.message, HTTPStatus.BAD_REQUEST)
        # TODO support CSV
        try:
            # TODO store mime_type when more than one is supported
            get_mime_type(acceptable_types=["application/json"], query_param=accept_type,
                          header=request.accept_mimetypes)
        except MimeTypeError as e:
            return make_response(e.message, HTTPStatus.NOT_ACCEPTABLE)
        try:
            return make_response((jsonify(current_app.data.data_frame(filter_, axis=Axis.OBS))), HTTPStatus.OK)
        except FilterError as e:
            return make_response(e.message, HTTPStatus.BAD_REQUEST)
        except ValueError as e:
            # JSON encoding failure, usually due to bad data
            warnings.warn(JSON_NaN_to_num_warning_msg)
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)

    @swagger.doc({
        "summary": "Get data (expression values) from the dataframe.",
        "tags": ["data"],
        "parameters": [
            {
                'name': 'filter',
                'description': 'Complex Filter',
                'in': 'body',
                'schema': FilterModel
            }
        ],
        "responses": {
            "200": {
                "description": "expression",
                "examples": {
                    "application/json": {
                        "var": [0, 20000],
                        "obs": [
                            [1, 39483, 3902, 203, 0, 0, 28]
                        ]
                    }
                }
            },
            "400": {
                "description": "Malformed filter"
            },
            "406": {
                "description": "Unacceptable MIME type"
            },
        }
    })
    def put(self):
        if not request.accept_mimetypes.best_match(["application/json", "text/csv"]):
            return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
        try:
            get_mime_type(acceptable_types=["application/json"], header=request.accept_mimetypes)
        except MimeTypeError as e:
            return make_response(e.message, HTTPStatus.NOT_ACCEPTABLE)
        try:
            return make_response((jsonify(current_app.data.data_frame(request.get_json()["filter"], axis=Axis.OBS))),
                                 HTTPStatus.OK)
        except FilterError as e:
            return make_response(e.message, HTTPStatus.BAD_REQUEST)
        except ValueError as e:
            # JSON encoding failure, usually due to bad data
            warnings.warn(JSON_NaN_to_num_warning_msg)
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


class DataVarAPI(Resource):
    @swagger.doc({
        "summary": "Get data (expression values) from the dataframe.",
        "tags": ["data"],
        "parameters": [
            {
                "in": "query",
                "name": "filter",
                "type": "string",
                "description": "axis:key:value"
            },
            {
                "in": "query",
                "name": "accept-type",
                "type": "string",
                "description": "MIME type"
            },
        ],
        "responses": {
            "200": {
                "description": "expression",
                "examples": {
                    "application/json": {
                        "obs": [0, 20000],
                        "var": [
                            [1, 39483, 3902, 203, 0, 0, 28]
                        ]
                    }
                }
            },
            "400": {
                "description": "Malformed filter"
            },
            "406": {
                "description": "Unacceptable MIME type"
            },
        }
    })
    def get(self):
        accept_type = request.args.get("accept-type", None)
        # request.args is immutable
        args = request.args.copy()
        args.pop("accept-type", None)
        try:
            filter_ = parse_filter(ImmutableMultiDict(args), current_app.data.schema['annotations'])
        except QueryStringError as e:
            return make_response(e.message, HTTPStatus.BAD_REQUEST)
        try:
            get_mime_type(acceptable_types=["application/json"], query_param=accept_type,
                          header=request.accept_mimetypes)
        except MimeTypeError as e:
            return make_response(e.message, HTTPStatus.NOT_ACCEPTABLE)
        try:
            return make_response((jsonify(current_app.data.data_frame(filter_, axis=Axis.VAR))), HTTPStatus.OK)
        except FilterError as e:
            return make_response(e.message, HTTPStatus.BAD_REQUEST)
        except ValueError as e:
            # JSON encoding failure, usually due to bad data
            warnings.warn(JSON_NaN_to_num_warning_msg)
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)

    @swagger.doc({
        "summary": "Get data (expression values) from the dataframe.",
        "tags": ["data"],
        "parameters": [
            {
                'name': 'filter',
                'description': 'Complex Filter',
                'in': 'body',
                'schema': FilterModel
            }
        ],
        "responses": {
            "200": {
                "description": "expression",
                "examples": {
                    "application/json": {
                        "obs": [0, 20000],
                        "var": [
                            [1, 39483, 3902, 203, 0, 0, 28]
                        ]
                    }
                }
            },
            "400": {
                "description": "Malformed filter"
            },
            "406": {
                "description": "Unacceptable MIME type"
            },
        }
    })
    def put(self):
        if not request.accept_mimetypes.best_match(["application/json", "text/csv"]):
            return make_response(f"Unsupported MIME type '{request.accept_mimetypes}'", HTTPStatus.NOT_ACCEPTABLE)
        # TODO support CSV
        try:
            get_mime_type(acceptable_types=["application/json"], header=request.accept_mimetypes)
        except MimeTypeError as e:
            return make_response(e.message, HTTPStatus.NOT_ACCEPTABLE)
        try:
            return make_response((jsonify(current_app.data.data_frame(request.get_json()["filter"], axis=Axis.VAR))),
                                 HTTPStatus.OK)
        except FilterError as e:
            return make_response(e.message, HTTPStatus.BAD_REQUEST)
        except ValueError as e:
            # JSON encoding failure, usually due to bad data
            warnings.warn(JSON_NaN_to_num_warning_msg)
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


class DiffExpObsAPI(Resource):
    @swagger.doc({
        "summary": "Generate differential expression (DE) statistics for two specified subsets of data, "
                   "as indicated by the two provided observation complex filters",
        "tags": ["diffexp"],
        # TODO sort out params
        # "parameters": [
        #     # {
        #     #     "in": "body",
        #     #     "name": "mode",
        #     #     "type": "string",
        #     #     "required": True,
        #     #     "description": "topN or varFilter"
        #     # },
        #     {
        #         "in": "query",
        #         "name": "count",
        #         "type": "int32",
        #         "description": "TopN mode: how many vars to return"
        #     },
        #     {
        #         "in": "body",
        #         "name": "varFilter",
        #         "schema": FilterModel,
        #         "description": "varFilter: Complex filter, only var for which vars to return"
        #     },
        #     {
        #         "in": "body",
        #         "name": "set1",
        #         "schema": FilterModel,
        #         "required": True,
        #         "description": "Complex filter, only obs - observations in set1"
        #     },
        #     {
        #         "in": "body",
        #         "name": "set2",
        #         "schema": FilterModel,
        #         "description": "Complex filter, only obs - observations in set2. If not included, inverse of set1."
        #     },
        # ],
        "responses": {
            "200": {
                "description": "Statistics are encoded as an array of arrays, with fields ordered as: "
                               "varIndex, logfoldchange,  pVal, pValAdj",
                "examples": {
                    "application/json": [
                        [328, -2.569489, 2.655706e-63, 3.642036e-57],
                        [1250, -2.569489, 2.655706e-63, 3.642036e-57],
                    ]
                }
            },
            "400": {
                "description": "malformed filter"
            },
            "403": {
                "description": "non-interactive request"
            },
            "501": {
                "description": "diffexp is not implemented"
            }
        }
    })
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
            diffexp = current_app.data.diffexp_topN(set1_filter, set2_filter, count,
                                                    current_app.data.features["diffexp"]["interactiveLimit"])
        except (ValueError, FilterError) as e:
            return make_response(e.message, HTTPStatus.BAD_REQUEST)
        except InteractiveError:
            return make_response("Non-interactive request", HTTPStatus.FORBIDDEN)
        try:
            return make_response(jsonify(diffexp), HTTPStatus.OK)
        except ValueError as e:
            # JSON encoding failure, usually due to bad data
            warnings.warn(JSON_NaN_to_num_warning_msg)
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)


class LayoutObsAPI(Resource):
    @swagger.doc({
        "summary": "Get the default layout for all observations.",
        "tags": ["layout"],
        "parameters": [],
        "responses": {
            "200": {
                "description": "layout",
                "examples": {
                    "application/json": {
                        "layout": {
                            "ndims": 2,
                            "coordinates": [
                                [0, 0.284483, 0.983744],
                                [1, 0.038844, 0.739444]
                            ]
                        }
                    }
                }
            },
            "400": {
                "description": "Data preparation error"
            }
        }
    })
    def get(self):
        try:
            layout = current_app.data.layout({})
        except PrepareError as e:
            return make_response(e.message, HTTPStatus.INTERNAL_SERVER_ERROR)
        try:
            return make_response((jsonify({"layout": layout})), HTTPStatus.OK)
        except ValueError as e:
            # JSON encoding failure, usually due to bad data
            warnings.warn(JSON_NaN_to_num_warning_msg)
            return make_response(str(e), HTTPStatus.INTERNAL_SERVER_ERROR)

    # @swagger.doc({
    #     "summary": "Observation layout for filtered subset.",
    #     "tags": ["layout"],
    #     "parameters": [
    #         {
    #             "name": "filter",
    #             "description": "Complex Filter",
    #             "in": "body",
    #             "schema": FilterModel
    #         }
    #     ],
    #     "responses": {
    #         "200": {
    #             "description": "layout",
    #             "examples": {
    #                 "application/json": {
    #                     "layout": {
    #                         "ndims": 2,
    #                         "coordinates": [
    #                             [0, 0.284483, 0.983744],
    #                             [1, 0.038844, 0.739444]
    #                         ]
    #                     }
    #                 }
    #             }
    #         },
    #         "400": {
    #             "description": "Malformed filter"
    #         },
    #         "403": {
    #             "description": "Non-interactive request"
    #         },
    #     }
    # })
    # def put(self):
    #     try:
    #         filter = request.get_json()["filter"]
    #         interactive_limit = current_app.data.features["layout"]["obs"]["interactiveLimit"]
    #         layout = current_app.data.layout(filter, interactive_limit=interactive_limit)
    #         return make_response(jsonify({"layout": layout}), HTTPStatus.OK)
    #     except FilterError as e:
    #         return make_response(e.message, HTTPStatus.BAD_REQUEST)
    #     except InteractiveError:
    #         return make_response("Non-interactive request", HTTPStatus.FORBIDDEN)


def get_api_resources():
    bp = Blueprint("api", __name__, url_prefix="/api/v0.2")
    api = Api(bp, add_api_spec_resource=False)
    # Initialization routes
    api.add_resource(SchemaAPI, "/schema")
    api.add_resource(ConfigAPI, "/config")
    # Data routes
    api.add_resource(AnnotationsObsAPI, "/annotations/obs")
    api.add_resource(AnnotationsVarAPI, "/annotations/var")
    api.add_resource(DataObsAPI, "/data/obs")
    api.add_resource(DataVarAPI, "/data/var")
    # Computation routes
    api.add_resource(DiffExpObsAPI, "/diffexp/obs")
    api.add_resource(LayoutObsAPI, "/layout/obs")
    return api
