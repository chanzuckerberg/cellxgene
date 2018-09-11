from flask_restful_swagger_2 import Schema


class AnnotationModel(Schema):
    type = "object"
    description = "Filter by annotation key: value"
    properties = {
        "name": {
            "type": "string"
        },
        # TODO update to OpenAPI v3.0 when a library is available that supports it
        # Unfortunately 2.0 doesn't have a way to have a schema that accepts multiple types
        # Overloading the type key with a list seems to work ok and makes it to the page
        "values": {
            "type": "array",
            "items": {
                "type": ["float32", "string", "int32", "bool"]
            }
        },
        "min": {
            "type": ["int32", "float32"],
        },
        "max": {
            "type": ["int32", "float32"],
        }
    }
    required = ["name"]


class IndexModel(Schema):
    type = "object"
    description = "Filter by index of observation/variable ex. [0, 5, 15]"
    properties = {
        "index": {
            "type": "array",
            "items": {
                "format": "int32",
                "type": "integer"
            }

        }
    }


class AxisModel(Schema):
    type = "object"
    description = "Axis of data -- obs or var"
    properties = {
        "index": IndexModel,
        "annotation_value": AnnotationModel.array()
    }


class FilterModel(Schema):
    type = "object"
    description = "Complex filter"
    properties = {
        "filter": {
            "type": "object",
            "properties": {
                "obs": AxisModel,
                "var": AxisModel
            }
        }
    }
