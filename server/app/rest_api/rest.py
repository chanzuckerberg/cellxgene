import pkg_resources

from flask import (
    Blueprint, current_app, jsonify, make_response, request
)
from flask_restful_swagger_2 import Api, swagger, Resource

from server.app.util.models import FilterModel
from server.app.util.constants import Axis, DiffExpMode


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
        return make_response(jsonify({"schema": current_app.data.schema}), 200)


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
                }
            }
        }
        return make_response(jsonify(config), 200)


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
            }
        }
    })
    def get(self):
        return make_response((jsonify({"layout": current_app.data.layout(current_app.data.data)})))

    @swagger.doc({
        "summary": "Observation layout for filtered subset.",
        "tags": ["layout"],
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
            }
        }
    })
    def put(self):
        df = current_app.data.filter_dataframe(request.get_json()["filter"])
        return make_response((jsonify({"layout": current_app.data.layout(df)})))


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
                            'tissue_type', 'sex', 'num_reads', 'clusters'
                        ],
                        "data": [
                            [0, 'lung', 'F', 39844, 99],
                            [1, 'heart', 'M', 83, 1],
                            [49, 'spleen', None, 2, "unknown cluster"],

                        ]
                    }

                }
            }
        }
    })
    def get(self):
        fields = request.args.getlist("annotation-name", None)
        try:
            annotation_response = current_app.data.annotation(current_app.data.data, fields)
        except KeyError:
            return make_response(f"Error bad key in {fields}", 404)
        else:
            return make_response(jsonify(annotation_response))

    @swagger.doc({
        "summary": "Fetch annotations (metadata) for filtered subset.",
        "tags": ["annotations"],
        "parameters": [
            {
                "in": "query",
                "name": "annotation-name",
                "type": "string",
                "description": "list of 1 or more annotation names"
            },
            {
                'name': 'filter',
                'description': 'Complex Filter',
                'in': 'body',
                'schema': FilterModel
            }
        ],
        "responses": {
            "200": {
                "description": "annotations",
                "examples": {
                    "application/json": {
                        "names": [
                            'tissue_type', 'sex', 'num_reads', 'clusters'
                        ],
                        "data": [
                            [0, 'lung', 'F', 39844, 99],
                            [1, 'heart', 'M', 83, 1],
                            [49, 'spleen', None, 2, "unknown cluster"],

                        ]
                    }

                }
            }
        }
    })
    def put(self):
        fields = request.args.getlist("annotation-name", None)
        df = current_app.data.filter_dataframe(request.get_json()["filter"])
        try:
            annotation_response = current_app.data.annotation(df, fields)
        except KeyError:
            return make_response(f"Error bad key in {fields}", 404)
        return make_response(jsonify(annotation_response))


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
                               "varIndex, avgDiff,  pVal, pValAdj, set1AvgExp, set2AvgExp",
                "examples": {
                    "application/json": [
                        [328, -2.569489, 2.655706e-63, 3.642036e-57, 383.393, 583.9],
                        [1250, -2.569489, 2.655706e-63, 3.642036e-57, 383.393, 583.9],
                    ]
                }
            }
        }
    })
    def post(self):
        args = request.get_json()
        mode = args["mode"]
        # Validate filters
        if "varFilter" in args:
            if mode != DiffExpMode.VAR_FILTER:
                return make_response(f"Var filter included but mode requested is {mode} not varFilter")
            if Axis.OBS in args["varFilter"]["filter"]:
                return make_response("Obs filter not allowed in varFilter", 400)
        # set1 filter only obs
        if Axis.VAR in args["set1"]["filter"]:
            return make_response("Var filter not allowed for set1", 400)
        # set2
        if "set2" in args and Axis.VAR in args["set2"]["filter"]:
            return make_response("Var filter not allowed for set2", 400)
        if mode == DiffExpMode.TOP_N and "set2" not in args:
            return make_response("Set2 as inverse of set1 is not implemented", 501)
        set1Filter = args["set1"]["filter"]
        set2Filter = args.get("set2", {"filter": {}})["filter"]
        if "varFilter" in args:
            set1Filter[Axis.VAR] = args["varFilter"]["filter"][Axis.VAR]
            set2Filter[Axis.VAR] = args["varFilter"]["filter"][Axis.VAR]
        df1 = current_app.data.filter_dataframe(set1Filter)
        # TODO inverse?
        df2 = current_app.data.filter_dataframe(set2Filter)
        # TODO exceeds size limit
        # mode
        count = args.get("count", None)
        try:
            diffexp = current_app.data.diffexp(df1, df2, count)
        except ValueError as ve:
            return make_response(ve.message, 400)
        return make_response(jsonify(diffexp))


def get_api_resources():
    bp = Blueprint("api", __name__, url_prefix="/api/v0.2")
    api = Api(bp, add_api_spec_resource=False)
    api.add_resource(SchemaAPI, "/schema")
    api.add_resource(ConfigAPI, "/config")
    api.add_resource(LayoutObsAPI, "/layout/obs")
    api.add_resource(AnnotationsObsAPI, "/annotations/obs")
    return api
