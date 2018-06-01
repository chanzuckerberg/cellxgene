class QueryStringError(Exception):
    pass

def _convert_variable(datatype, variable):
    """
    Convert variable to number (float/int)
    Used for dataset metadata and for query string
    :param datatype: type to convert to
    :param variable: value of variable
    :return: converted variable
    :raises: ValueError
    """
    try:
        if variable and datatype == "int":
            variable = int(variable)
        elif variable and datatype == "float":
            variable = float(variable)
        return variable
    except ValueError:
        raise

def parse_filter(filter, schema):
    """
    {key: variable_type
    value_type
    query

    :param filter:
    :param schema:
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
            query[key]["query"] = _convert_variable(query[key]["value_type"], value)
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
