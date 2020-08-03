import warnings

import numpy as np
import pandas as pd


def series_to_schema(array: pd.Series):
    try:
        return dtype_to_schema(array.dtype)
    except TypeError:
        dtype = array.dtype
        data_kind = dtype.kind
        schema = {}
        if can_cast_to_float32(array):
            schema["type"] = "float32"
        elif can_cast_to_int32(array):
            schema["type"] = "int32"
        elif data_kind == "O" and dtype == "object":
            schema["type"] = "string"
        else:
            raise TypeError(f"Annotations of type {dtype} are unsupported.")
        return schema


def dtype_to_schema(dtype):
    schema = {}
    if dtype == np.float32:
        schema["type"] = "float32"
    elif dtype == np.int32:
        schema["type"] = "int32"
    elif dtype == np.bool_:
        schema["type"] = "boolean"
    elif dtype == np.str:
        schema["type"] = "string"
    elif dtype == "category":
        schema["type"] = "categorical"
        schema["categories"] = dtype.categories.tolist()
    else:
        raise TypeError(f"Annotations of type {dtype} are unsupported.")
    return schema


def can_cast_to_float32(array):
    if array.dtype.kind == "f":
        if not np.can_cast(array.dtype, np.float32):
            warnings.warn(f"Annotation {array.name} will be converted to 32 bit float and may lose precision.")
        return True
    return False


def can_cast_to_int32(array):
    if array.dtype.kind in ["i", "u"]:
        if np.can_cast(array.dtype, np.int32):
            return True
        ii32 = np.iinfo(np.int32)
        if array.min() >= ii32.min and array.max() <= ii32.max:
            return True
    return False
