class QueryStringError(Exception):
    pass


def _convert_variable(datatype, variable):
    """
    Convert variable to number (float/int)
    Used for dataset metadata and for query string
    :param datatype: type to convert to
    :param variable (string or None): value of variable
    :return: converted variable
    :raises: ValueError
    """
    try:
        if variable is None:
            return variable
        if datatype == "int":
            variable = int(variable)
        elif datatype == "float":
            variable = float(variable)
        return variable
    except ValueError:
        raise


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
    :return:
    """
    query = {}
    for key in filter:
        value = filter.getlist(key)
        if key not in schema:
            raise QueryStringError("Error: key {} not in metadata schema".format(key))
        query[key] = {
            "variable_type": schema[key]["variabletype"],
            "value_type": schema[key]["type"]
        }
        if query[key]["variable_type"] == "categorical":
            query[key]["query"] = [_convert_variable(query[key]["value_type"], v) for v in value]
        elif query[key]["variable_type"] == "continuous":
            value = value[0]
            try:
                min, max = value.split(",")
            except ValueError:
                raise QueryStringError("Error: min,max format required for range for key {}, got {}".format(key, value))
            if min == "*":
                min = None
            if max == "*":
                max = None
            try:
                query[key]["query"] = {
                    "min": _convert_variable(query[key]["value_type"], min),
                    "max": _convert_variable(query[key]["value_type"], max)
                }
            except ValueError:
                raise QueryStringError(
                    "Error: expected type {} for key {}, got {}".format(query[key]["type"], key, value)
                )
    return query
