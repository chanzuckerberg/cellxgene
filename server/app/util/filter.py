from numpy import float32, int32

from server.app.util.constants import Axis


class QueryStringError(Exception):
    pass


def _convert_variable(datatype, variable):
    """
    Convert variable to number (float/int)
    Used for dataset metadata and for query string
    :param datatype: type to convert to
    :param variable (string or None): value of variable
    :return: converted variable
    :raises: AssertionError, ValueError
    """
    assert datatype in ["boolean", "categorical", "float32", "int32", "string"]
    if variable is None:
        return variable
    try:
        if datatype == "int32":
            variable = int32(variable)
        elif datatype == "float32":
            variable = float32(variable)
        elif datatype == "boolean":
            if variable == "false":
                variable = False
            elif variable == "true":
                variable = True
            else:
                raise ValueError
    except ValueError:
        raise
    else:
        return variable


def parse_filter(filter, schema):
    """
    The filter comes in as arguments from a GET/POST request
    For categorical metadata keys filter based on key=value
    For continuous metadata keys filter by key=min,max
    Either value can be replaced by a * To have only a minimum value key=min, To have only a maximum value key=*,max

    They combine via AND so a cell's metadata would have to match every filter

    The results is a matrix with the cells the pass the filter and at this point all the genes
    :param filter: flask's request.args
    :param schema: dictionary schema
    :raises QueryStringError
    :return:
    """
    query = {}

    for key in filter:
        axis, annotation = key.split(":", 1)
        if axis not in [ax for ax in Axis]:
            raise QueryStringError(f"Error: bad axis in query string {axis}")
        if axis not in query:
            query[axis] = {}
        if annotation == "index":
            query[axis]["index"] = filter.getlist(key)
        else:
            if "annotation_value" not in query[axis]:
                query[axis]["annotation_value"] = []
            ann_filter = {"name": annotation}
            for ann in schema[axis]:
                if ann["name"] == annotation:
                    dtype = ann["type"]
                    break
            else:
                raise QueryStringError(f"Error: {annotation} not a valid annotation type")
            if dtype in ["string", "categorical", "boolean"]:
                ann_filter["value"] = [_convert_variable(dtype, i) for i in filter.getlist(key)]
            else:
                value = filter.get(key)
                try:
                    min_, max_ = value.split(",")
                except ValueError:
                    raise QueryStringError(f"Error: min,max format required for range for {annotation}, got {value}")
                if min_ == "*":
                    min_ = None
                if max_ == "*":
                    max_ = None
                try:
                    ann_filter["value"] = {
                        "min": _convert_variable(dtype, min_),
                        "max": _convert_variable(dtype, max_)
                    }
                except ValueError:
                    raise QueryStringError(
                        f"Error: expected type {dtype} for key {annotation}, got {value}"
                    )
            query[axis]["annotation_value"].append(ann_filter)
    return query
