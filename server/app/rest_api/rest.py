from flask import (
    Blueprint, current_app, jsonify, make_response
)
from flask_restful_swagger_2 import Api, swagger, Resource


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
                    "engine": f"ScanPy Engine version {current_app.data.__version__}",
                    "dataset": current_app.config["DATASET_TITLE"]
                }
            }
        }
        return make_response(jsonify(config), 200)


def get_api_resources():
    bp = Blueprint("api", __name__, url_prefix="/api/v0.2")
    api = Api(bp, add_api_spec_resource=False)
    api.add_resource(SchemaAPI, "/schema")
    api.add_resource(ConfigAPI, "/config")
    return api
